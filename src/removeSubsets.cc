// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 


/* jonathan myers this is an attempt at doing removal of subset tracklets in a
 * much faster way than the naive algorithm.
 *
 * someday, it might make sense to template-ize this...
*/


#include <algorithm> /* for a decent set intersection implementation...*/
#include <fstream>
#include <iostream>
#include <map>
#include <set>

#include "lsst/mops/removeSubsets.h"
#include "lsst/mops/Exceptions.h"




namespace lsst {
    namespace mops {

    void SubsetRemover::removeSubsetsPopulateOutputVector(const std::vector<Tracklet> *pairsVector, 
                                                          std::vector<Tracklet> &outVector) {
        /* build a map which maps each detection to each tracklet which uses it. */
        
        std::map<unsigned int, std::set<unsigned int> > detectionIndexToPairsVectorIndexSet;
        std::set<unsigned int>::iterator indicesIter;
        
        for (unsigned int curPairIndex = 0; curPairIndex < pairsVector->size(); curPairIndex++) {

            const Tracklet *curPair = &(*pairsVector)[curPairIndex];

            for (indicesIter = curPair->indices.begin(); 
                 indicesIter != curPair->indices.end();
                 indicesIter++) {
                /* indicesIter now points to an index into the detections file (i.e. detection ID) */
                detectionIndexToPairsVectorIndexSet[*indicesIter].insert(curPairIndex);
            }
        }
        
        std::vector<bool> pairHasBeenWritten(pairsVector->size(), false); 

        /* for each tracklet: see if any other tracklet uses all the same
         * detections as this one.  If it does, it is either a superset or an
         * identical tracklet. */

        for (unsigned int curPairIndex = 0; curPairIndex < pairsVector->size(); curPairIndex++)  {

            const Tracklet * curPair =  &(*pairsVector)[curPairIndex];
            /* iteratively calculate the intersect of the various sets of pairs
             * which use the detections in this pair.  Start with the first set
             * - the one associated with the first detection in this
             * tracklet.
             *
             *it is safe to assume that all tracklets have at least 1
             * detection, of course */
            std::set<unsigned int> indicesIntersect = detectionIndexToPairsVectorIndexSet[*(curPair->indices.begin())]; 

            std::set<unsigned int>::const_iterator indexIter;
            /* now go ahead and try iterating through the rest of the results... */
            for (indexIter = (curPair->indices.begin())++; 
                 indexIter != curPair->indices.end(); 
                 indexIter++) {
                std::set<unsigned int> tmpResults;
                const std::set<unsigned int> * detIDSet = &(detectionIndexToPairsVectorIndexSet[*indexIter]);
                std::set_intersection(indicesIntersect.begin(), indicesIntersect.end(),
                                      detIDSet->begin(), detIDSet->end(), 
                                      std::insert_iterator<std::set <unsigned int> >(tmpResults, tmpResults.begin()));
                indicesIntersect = tmpResults;                                      
            }
            /* indicesIntersect now holds the set of all indices into
             * pairsVector s.t. the tracklet at that index holds every element in
             * the current tracklet.*/

            if (indicesIntersect.size() == 1) {
                /* we are the only tracklet with these detections. */
                pairHasBeenWritten[curPairIndex] = true;
                outVector.push_back(*curPair);                
            }
            else if (indicesIntersect.size() > 1) {
                /* we got other tracklet(s) which may be supersets or identical.  Check for supersets.*/
                bool writeThisTracklet = true;
                unsigned int mySize = curPair->indices.size();
                std::set<unsigned int>::iterator otherTrackletIter;

                for (otherTrackletIter = indicesIntersect.begin(); 
                     otherTrackletIter != indicesIntersect.end();
                     otherTrackletIter++) {

                    if ((*pairsVector)[*otherTrackletIter].indices.size() > mySize) {
                        writeThisTracklet = false;
                    }

                    else if (((*pairsVector)[*otherTrackletIter].indices.size() == mySize) && 
                             (pairHasBeenWritten[*otherTrackletIter] == true)) {
                        writeThisTracklet = false;
                    }

                }
                if (writeThisTracklet == true) {
                    /* we had identical tracklets, but they haven't been written
                     * out - so write ourselves, and mark ourselves as
                     * written */                    
                    pairHasBeenWritten[curPairIndex] = true;
                    outVector.push_back(*curPair);                
                }

            }
            else {
                throw LSST_EXCEPT(ProgrammerErrorException, "Programming error: impossible case in removeSubsets.cc");
            }
            
        }
    }



    



    void putLongestPerDetInOutputVector(const std::vector<Tracklet> *pairsVector, 
                                        std::vector<Tracklet> &outputVector) {

        std::map<unsigned int, std::set<unsigned int> > detectionIndexToPairsVectorIndexSet;
        std::set<unsigned int>::iterator indicesIter;

        /* make a mapping from each detID to the tracklet IDs which contain that det */

        for (unsigned int curPairIndex = 0; curPairIndex < pairsVector->size(); curPairIndex++) {
            
            const Tracklet *curPair = &(*pairsVector)[curPairIndex];

            for (indicesIter = curPair->indices.begin(); 
                 indicesIter != curPair->indices.end();
                 indicesIter++) {
                /* indicesIter now points to an index into the detections file (i.e. detection ID) */
                detectionIndexToPairsVectorIndexSet[*indicesIter].insert(curPairIndex);
            }
        }
        /* now for each det, find the longest associated tracklet.  Add that one
         * to outputVector. (if not already added). */

        std::vector<bool> writeThisIndex(pairsVector->size(), false);

        /* for each detection... */
        std::map<unsigned int, std::set<unsigned int> >::iterator detIter;
        for (detIter = detectionIndexToPairsVectorIndexSet.begin();
             detIter != detectionIndexToPairsVectorIndexSet.end();
             detIter++) {
            unsigned int maxLen = 0;
            unsigned int longestTrackletIndex = 0;
            /* for each tracklet associated with that detection... */
            for (indicesIter = detIter->second.begin(); 
                 indicesIter != detIter->second.end();
                 indicesIter++) {
                /*indicesiter points at an int which is an index into pairsVector. */
                unsigned int curTrackletLen = (*pairsVector)[*indicesIter].indices.size();
                if ( curTrackletLen > maxLen) {
                    maxLen = curTrackletLen;
                    longestTrackletIndex = *indicesIter;
                }
            }
            /* make note that this tracklet was chosen */
            if (maxLen > 0) {
                writeThisIndex[longestTrackletIndex] = true;
            }
        }
        /* now write out chosen tracklet to output vector */
        for (unsigned int i = 0; i < writeThisIndex.size(); i++) {
            if (writeThisIndex[i] == true) {
                outputVector.push_back((*pairsVector)[i]);
            }
        }
    }








    }} // close lsst::mops



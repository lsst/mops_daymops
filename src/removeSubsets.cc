// -*- LSST-C++ -*-


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
#include <time.h>

#include "lsst/mops/removeSubsets.h"
#include "lsst/mops/Exceptions.h"




namespace lsst {
    namespace mops {

std::string curTime()
{
    time_t timer;
    timer=time(NULL);
    return std::string(asctime(localtime(&timer)));
}


bool trackContains(unsigned int ind, const Tracklet &t) 
{
    return (t.indices.find(ind) != t.indices.end());
}




void getSupersetOrIdenticalTracks(const Tracklet* curTrack, 
                                  std::map<unsigned int, std::set<unsigned int> > & reverseMap,
                                  bool shortCircuit,
                                  bool sortBeforeIntersect,
                                  std::set<unsigned int> & results)
{

        
    std::set<unsigned int>::const_iterator indexIter;
    bool quitNow = false;
    if (!sortBeforeIntersect) 
    {
        // initialize results set
        results = reverseMap[*(curTrack->indices.begin())]; 

        // if not sorting, read elements from reverse-map as needed;
        // if we short-circuit we can avoid some memory accesses.
        for (indexIter = (curTrack->indices.begin())++; 
             indexIter != curTrack->indices.end() && (quitNow == false);
             indexIter++) {
            std::set<unsigned int> tmpResults;
            const std::set<unsigned int> * detIDSet = &(reverseMap[*indexIter]);
            std::set_intersection(results.begin(), results.end(),
                                  detIDSet->begin(), detIDSet->end(), 
                                  std::insert_iterator<std::set <unsigned int> >(tmpResults, tmpResults.begin()));
            results = tmpResults;
            
            if (shortCircuit && (results.size() == 1)) {
                // we are not a subset track(let); end this for loop
                quitNow = true;
            }
            
        }
    }


    else {

        /* we're sorting the sets read from reverseMap by size before
         * doing intersection; so first we need to read all the
         * elements from the reverse map, sort them by size and THEN
         * do the intersection */
        std::multimap<unsigned int, const std::set<unsigned int>* > reverseMapEntries;
        for (indexIter = (curTrack->indices.begin()); 
             indexIter != curTrack->indices.end();
             indexIter++) {
            const std::set<unsigned int>* detIDSet = &(reverseMap[*indexIter]);
            std::pair<unsigned int, const std::set<unsigned int> * > toInsert(detIDSet->size(), detIDSet);
            reverseMapEntries.insert(toInsert);
        }

        /* std::multimap sorts automatically so we'll be using smallest
           element first */
        std::multimap<unsigned int, const std::set<unsigned int>* >::const_iterator entryIter;
        // initialize results set
        results = *(reverseMapEntries.begin()->second);

        for (entryIter = reverseMapEntries.begin()++; 
             entryIter != reverseMapEntries.end() && (quitNow == false);
             entryIter++) {

            std::set<unsigned int> tmpResults;
            std::set_intersection(results.begin(), results.end(),
                                  entryIter->second->begin(), entryIter->second->end(), 
                                  std::insert_iterator<std::set <unsigned int> >(tmpResults, tmpResults.begin()));
            results = tmpResults;
            
            if (shortCircuit && (results.size() == 1)) {
                // we are not a subset track(let); end this for loop
                quitNow = true;
            }
            
        }
           
    }


}







void SubsetRemover::removeSubsetsPopulateOutputVector(const std::vector<Tracklet> *tracksVector, 
                                                      std::vector<Tracklet> &outVector,
                                                      bool shortCircuit, 
                                                      bool sortBeforeIntersect) {
    /* build a map which maps each detection to each tracklet which uses it. */
    
    std::map<unsigned int, std::set<unsigned int> > reverseMap;
    std::set<unsigned int>::iterator indicesIter;
    
    std::cout << "Building detection-to-track map, starting at " << curTime() << std::endl;
    for (unsigned int curTrackIndex = 0; curTrackIndex < tracksVector->size(); curTrackIndex++) {
        
        const Tracklet *curTrack = &(*tracksVector)[curTrackIndex];
        
        for (indicesIter = curTrack->indices.begin(); 
             indicesIter != curTrack->indices.end();
             indicesIter++) {
            /* indicesIter now points to an index into the detections file (i.e. detection ID) */
            reverseMap[*indicesIter].insert(curTrackIndex);
        }
    }
        
    std::cout << "Finished detection-to-track map, filtering tracks starting at " << curTime() << std::endl;

    std::vector<bool> trackHasBeenWritten(tracksVector->size(), false); 

    /* for each tracklet: see if any other tracklet uses all the same
     * detections as this one.  If it does, it is either a superset or an
     * identical tracklet. */

    for (unsigned int curTrackIndex = 0; curTrackIndex < tracksVector->size(); curTrackIndex++)  {
 
        if ((curTrackIndex > 0) && (curTrackIndex % 10000 == 0)) {
            std::cout << "Processing element " << curTrackIndex << " of " << tracksVector->size() << " (" 
                      << 100. * (curTrackIndex*1.0) / (tracksVector->size() * 1.0) << "%)" << std::endl;
        }

        const Tracklet * curTrack =  &(*tracksVector)[curTrackIndex];
        /* iteratively calculate the intersect of the various sets of tracks
         * which use the detections in this track.  Start with the first set
         * - the one associated with the first detection in this
         * tracklet.
         *
         *it is safe to assume that all tracklets have at least 1
         * detection, of course */

        std::set<unsigned int> indicesIntersect;  
        getSupersetOrIdenticalTracks(curTrack, 
                                     reverseMap,
                                     shortCircuit,
                                     sortBeforeIntersect, 
                                     indicesIntersect);


        /* indicesIntersect now holds the set of all indices into
         * tracksVector s.t. the tracklet at that index holds every element in
         * the current tracklet.*/

        if (indicesIntersect.size() == 1) {
            /* we are the only tracklet with these detections. */
            trackHasBeenWritten[curTrackIndex] = true;
            outVector.push_back(*curTrack);                
        }
        else if (indicesIntersect.size() > 1) {
            /* we got other tracklet(s) which may be supersets or identical.  Check for supersets.*/
            bool writeThisTracklet = true;
            unsigned int mySize = curTrack->indices.size();
            std::set<unsigned int>::iterator otherTrackletIter;

            for (otherTrackletIter = indicesIntersect.begin(); 
                 otherTrackletIter != indicesIntersect.end();
                 otherTrackletIter++) {

                if ((*tracksVector)[*otherTrackletIter].indices.size() > mySize) {
                    writeThisTracklet = false;
                }

                else if (((*tracksVector)[*otherTrackletIter].indices.size() == mySize) && 
                         (trackHasBeenWritten[*otherTrackletIter] == true)) {
                    writeThisTracklet = false;
                }

            }
            if (writeThisTracklet == true) {
                /* we had identical tracklets, but they haven't been written
                 * out - so write ourselves, and mark ourselves as
                 * written */                    
                trackHasBeenWritten[curTrackIndex] = true;
                outVector.push_back(*curTrack);                
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



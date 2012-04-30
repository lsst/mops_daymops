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
#include <omp.h>

#include "lsst/mops/removeSubsetsOMP.h"
#include "lsst/mops/Exceptions.h"

#define uint unsigned int



namespace lsst {
    namespace mops {

std::string curTime()
{
    time_t timer;
    timer=time(NULL);
    return std::string(asctime(localtime(&timer)));
}


bool trackContains(unsigned int ind, const std::set<unsigned int> &t) 
{
    return (t.find(ind) != t.end());
}




void getSupersetOrIdenticalTracks(const std::set<unsigned int>* curTrack, 
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
        results = reverseMap[*(curTrack->begin())]; 

        // if not sorting, read elements from reverse-map as needed;
        // if we short-circuit we can avoid some memory accesses.
        for (indexIter = (curTrack->begin())++; 
             indexIter != curTrack->end() && (quitNow == false);
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
        for (indexIter = (curTrack->begin()); 
             indexIter != curTrack->end();
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







void SubsetRemover::removeSubsetsPopulateOutputVector(
    const std::vector<std::set<unsigned int> > *tracksVector, 
    std::vector<bool> &outVector,
    bool shortCircuit,
    bool sortBeforeIntersect)
{
    /* build a map which maps each detection to each tracklet which uses it. */
    
    std::map<unsigned int, std::set<unsigned int> > reverseMap;
    std::set<unsigned int>::iterator indicesIter;
    

    int nthreads, tid;
    
#pragma omp parallel private(nthreads, tid)
    {
      nthreads = omp_get_num_threads();
      tid = omp_get_thread_num();
      if(tid == 0) {
	std::cout << "Number of threads " << nthreads << std::endl;
      }
    }

    std::cout << "Building detection-to-track map, starting at " << curTime() << std::endl;
    for (unsigned int curTrackIndex = 0; curTrackIndex < tracksVector->size(); curTrackIndex++) {
        
        const std::set<unsigned int> *curTrack = &(*tracksVector)[curTrackIndex];
        
        for (indicesIter = curTrack->begin(); 
             indicesIter != curTrack->end();
             indicesIter++) {
            /* indicesIter now points to an index into the detections file (i.e. detection ID) */
            reverseMap[*indicesIter].insert(curTrackIndex);
        }
    }
        
    std::cout << "Finished detection-to-track map, filtering tracks starting at " << curTime() << std::endl;

    outVector.clear();
    outVector.resize(tracksVector->size(), false); 

    /* for each tracklet: see if any other tracklet uses all the same
     * detections as this one.  If it does, it is either a superset or an
     * identical tracklet. */

#pragma omp parallel for schedule(dynamic, 1000)
    for (unsigned int curTrackIndex = 0; curTrackIndex < tracksVector->size(); curTrackIndex++)  {
 
        if ((curTrackIndex > 0) && (curTrackIndex % 10000 == 0)) {
            std::cout << "Processing element " << curTrackIndex << " of " << tracksVector->size() << " (" 
                      << 100. * (curTrackIndex*1.0) / (tracksVector->size() * 1.0) << "%)" << std::endl;
        }

        const std::set<unsigned int> * curTrack =  &(*tracksVector)[curTrackIndex];
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
            outVector.at(curTrackIndex) = true;
        }
        else if (indicesIntersect.size() > 1) {
            /* we got other tracklet(s) which may be supersets or
             * identical.  Check for supersets.*/
            bool writeThisTrack = true;
            unsigned int mySize = curTrack->size();
            std::set<unsigned int>::iterator otherTrackIter;

            for (otherTrackIter = indicesIntersect.begin(); 
                 otherTrackIter != indicesIntersect.end();
                 otherTrackIter++) {

                if ((*tracksVector)[*otherTrackIter].size() > mySize) {
                    writeThisTrack = false;
                }

                /* check if we got a redundant track. if so, keep this
                 * track ONLY IF ITS INDEX IN THE INPUT VECTOR is
                 * minimal among the redundant tracks. */
                else if (((*tracksVector)[*otherTrackIter].size() == mySize) && 
                         (curTrackIndex > *otherTrackIter)) {
                    writeThisTrack = false;
                }

            }
            if (writeThisTrack == true) {
                /* we had identical tracklets, but they haven't been written
                 * out - so write ourselves, and mark ourselves as
                 * written */                    
                outVector[curTrackIndex] = true;
            }

        }
        else {
            throw LSST_EXCEPT(ProgrammerErrorException, "Programming error: impossible case in removeSubsets.cc");
        }
            
    }
}



  



    }} // close lsst::mops



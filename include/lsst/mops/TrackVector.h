// -*- LSST-C++ -*-
/* jonathan myers 
   4/28/11
*/


/*
  TrackVector is not TrackSet.

  TrackSet is used implemented with a real Set; it is useful for doing
  subset/intersection comparisons, useful for e.g. the linkTracklets
  unit tests where we want to know whether the found tracks were a
  superset of the expected true tracks.  

  TrackVector, however, avoids this overhead (and allows faster
  lookup) because it is implemented with a vector.  Also it doesn't
  allow any of the fancy behavior allowing it to act as a cache for
  output.  

  However, it *can* be populated from a file, and can write helpful
  per-track stats to a file (e.g. epoch, velocity, acceleration)
 */

#ifndef LSST_TRACKVECTOR_H
#define LSST_TRACKVECTOR_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include <vector>

#include "lsst/mops/Track.h"
#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/LinkageVector.h"

namespace lsst {
namespace mops {


class TrackVector : public LinkageVector{
public:

    TrackVector() {} ;

    /* note: inFile MUST be in the form of 0-indexed indices into
     * allDets, NOT diaSourceIDs! */
    void populateFromFile(std::string fileName, 
                          const std::vector<MopsDetection> &allDets);

    /* write just the diaIds for each track. */
    void writeTracksToFile(std::string outFileName); 
    
    /* use allDets to calculate the best-fit quadratic, determine
     * whether tracks are true, and write those stats alongside the
     * diaSources. */
    void writeTracksAndStatsToFile(std::string outFileName,
                                   const std::vector<MopsDetection> &allDets);

    ~TrackVector() {};

    void debugPrint();

    void push_back(const Track *newTrack) { contents.push_back(*newTrack);}

    unsigned int size() const { return contents.size(); }

    Track* at(unsigned int index) { return &contents.at(index);}
        

private:
    std::vector<Track> contents;
    void writeSet(std::ofstream * outFile, const std::set<unsigned int> &s);

};



}} // close namespace lsst::mops

#endif

// -*- LSST-C++ -*-
/* jonathan myers 
   8/10/09
*/

#include <set>


#ifndef LSST_TRACK_H
#define LSST_TRACK_H


class Track {
public:
  // the tracklets which were used to build this track, if any
  std::set<unsigned int> componentTrackletIndices;
  /* the set of detections associated with this track.  note that this is not
     merely the union of the detections associated with the tracklets
     associated with this track; if the track contains two tracklets which
     are in 'conflict' (contain multiple detections from the same image time)
     then only one of the possible detections is chosen.        
  */
  std::set<unsigned int> componentDetectionIndices;
};




#endif

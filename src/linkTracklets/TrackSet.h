// -*- LSST-C++ -*-
/* jonathan myers 
   4/08/10
*/

#ifndef LSST_TRACKSET_H
#define LSST_TRACKSET_H

#include <set>
#include "Track.h"

class TrackSet {
public:
     std::set<Track> componentTracks;
     void insert(const Track &newTrack) {
	  componentTracks.insert(newTrack);
     };

     bool operator==(const TrackSet other) const {
	  return (componentTracks == other.componentTracks);
     };

    bool operator!=(const TrackSet &other) const {
        return ! (*this == other);
    }

};





#endif

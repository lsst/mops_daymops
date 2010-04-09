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
    unsigned int size() const {
        return componentTracks.size();
    };
    
    bool isSubsetOf(const TrackSet &other) const {
        if (other.size() < this->size()) {
            return false;
        }
        else {
            // create a copy TrackSet
            std::set<Track>::const_iterator tIter;
            for (tIter = componentTracks.begin();
                 tIter != componentTracks.end();
                 tIter++) {
                if (other.componentTracks.find(*tIter) 
                    == other.componentTracks.end()) {
                    return false;
                }
            }
            return true;
        }
    };

     bool operator==(const TrackSet other) const {
	  return (componentTracks == other.componentTracks);
     };

    bool operator!=(const TrackSet &other) const {
        return ! (*this == other);
    };

};





#endif

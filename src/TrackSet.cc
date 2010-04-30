// -*- LSST-C++ -*-
/* jonathan myers 
   4/08/10
*/


#include "lsst/mops/TrackSet.h"

namespace lsst { 
namespace mops {


void TrackSet::insert(const Track &newTrack) {
    componentTracks.insert(newTrack);
}



unsigned int TrackSet::size() const {
    return componentTracks.size();
}



bool TrackSet::isSubsetOf(const TrackSet &other) const {
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
}



bool TrackSet::operator==(const TrackSet other) const {
    return (componentTracks == other.componentTracks);
}



bool TrackSet::operator!=(const TrackSet &other) const {
    return ! (*this == other);
}

}} // close namespace lsst::mops

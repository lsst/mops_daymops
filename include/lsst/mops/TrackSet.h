// -*- LSST-C++ -*-
/* jonathan myers 
   4/08/10
*/

#ifndef LSST_TRACKSET_H
#define LSST_TRACKSET_H

#include <set>
#include "Track.h"


namespace lsst {
namespace mops {


class TrackSet {
public:
    std::set<Track> componentTracks;

    void insert(const Track &newTrack);

    unsigned int size() const;
    
    bool isSubsetOf(const TrackSet &other) const;

     bool operator==(const TrackSet other) const;

    bool operator!=(const TrackSet &other) const;

};



}} // close namespace lsst::mops

#endif

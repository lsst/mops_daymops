// -*- LSST-C++ -*-
/* jonathan myers 
   4/08/10
*/

#ifndef LSST_TRACKLETVECTOR_H
#define LSST_TRACKLETVECTOR_H

#include <vector>

#include "Tracklet.h"


namespace lsst {
namespace mops {

class TrackletVector {
public:
    std::vector<Tracklet> componentTracklets;

    void push_back(const Tracklet &newTracklet);

    Tracklet at(unsigned int) const;

    unsigned int size() const;
    
    /*
      this will return true iff this the other vector contains every tracklet
      that we do.  ordering is ignored.  

      to make the implementation FAST, we will create a COPY of the other
      trackletVector's tracklets; be wary of memory usage.
     */
    bool isSubsetOf(const TrackletVector &other) const;

    /*
     * see warnings in Tracklet.h.  If two tracklet sets were generated from
     * different detection vectors, even if those vectors contain the same
     * detections in a different ordering, the results of these comparison
     * functions are NOT meaningful!
     */

    /*
     * two vectors are == if and only if they contain the same tracklets;
     * ordering is NOT considered.  For fast performance, all Tracklets in other
     * will be copied, then deleted, then all tracklets in this TrackletVector
     * will be copied, then deleted.
     */
    bool operator==(const TrackletVector other) const;

    bool operator!=(const TrackletVector &other) const;

};

}} // close namespace lsst::mops



#endif

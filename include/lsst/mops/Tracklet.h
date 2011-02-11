// -*- LSST-C++ -*-


/* jmyers 8/18/08
 */



#ifndef LSST_TRACKLET_H
#define LSST_TRACKLET_H
#include <set>
#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/Exceptions.h"

namespace lsst {
namespace mops {

class Tracklet {
public: 
    Tracklet();
    Tracklet(std::set <unsigned int> startIndices);

    /* a tracklet is just a collection of detections.  To conserve memory, we
       don't hold copies of those detections here, we hold references to them.
       In all my code, I represent a detection with its index from the
       Detections vector I'm working with.  Note that these indices ARE NOT
       NECESSARILY the same as the ID's of the Detections themselves.
       
       It is TOTALLY UP TO THE USER of this class to make sure that they
       associate Tracklets with the correct Detection vector.
     */

    std::set<unsigned int> indices;
    bool isCollapsed;    
    // these fields are used only by linkTracklets. linkTracklets
    // is responsible for setting them before reading.
    double velocityRA;
    double velocityDec;
    
    /* return start time. does NOT assume indices are assigned
     * chronologically. ASSUMES that you send in the same vector of detections
     * used to create this tracklet! */
    double getStartTime(std::vector<MopsDetection> dets);

    /*
     * The following operators consider only the indices set.  Note that if you
     * compare two tracklets with indices into different detection vectors, the
     * results are MEANINGLESS! It is up to the user to ensure that
     * when comparing tracklets t1 and t2 that if t1.indices contains id X, then
     * X in t2.indices has the same meaning.
     */

    Tracklet & operator= (const Tracklet &other) {
        indices = other.indices;
        return *this;
    }

    bool operator==(const Tracklet &other) const {
        bool toRet = indices == other.indices;
        return toRet ;
    }

    bool operator!=(const Tracklet &other) const {
        return ! (*this == other);
    }

    bool operator<(const Tracklet &other) const {
        bool toRet = indices < other.indices;
        return toRet;
        
    }


};


} } //close namespace lsst::mops

#endif

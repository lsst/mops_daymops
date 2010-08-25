// -*- LSST-C++ -*-


/* jmyers 8/18/08
 *
 * previously part of collapseTracklets, this is our incredibly simple Tracklet
 * class (could probably be made a struct.).  indices is just a set of
 * detectionIDs.  Supply your own detection IDs.
 *
 */



#ifndef LSST_TRACKLET_H
#define LSST_TRACKLET_H
#include <set>

namespace lsst {
namespace mops {

class Tracklet {
public: 
    Tracklet() { isCollapsed = false; velocityRA = 0; velocityDec = 0;}
    Tracklet(std::set <unsigned int> startIndices) { isCollapsed = false; indices = startIndices;}

    /* a tracklet is just a collection of detections.  To conserve memory, we
       don't hold those detections in-mem, we hold references to them.  In all
       my code, I represent a detection with its index from the Detections
       vector I'm working with.  This way we can easily index Note that these
       ID's ARE NOT NECESSARILY the same as the ID's of the Detections
       themselves.  
       
       This is inelegant, but otherwise we'd either store a lot more data than
       we want (hurting performance) or we'd have to store a reference to the
       generating detection vector here - which is NOT what we want to do in
       certain circumstances, such as subset removal, in which we need not ever
       look at the detections themselves, just the indices sets on the various
       tracklets.
     */

    std::set<unsigned int> indices;
    bool isCollapsed;    
    // these fields are used only by linkTracklets. linkTracklets
    // is responsible for setting them before reading.
    double velocityRA;
    double velocityDec;


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

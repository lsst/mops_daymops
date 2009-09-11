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

class Tracklet {
public: 
    Tracklet() { isCollapsed = false; velocityRA = 0; velocityDec = 0;}
    Tracklet(std::set <unsigned int> startIndices) { isCollapsed = false; indices = startIndices;}
    std::set<unsigned int> indices;
    bool isCollapsed;    
    // these fields are used only by linkTracklets. linkTracklets
    // is responsible for setting them before reading.
    double velocityRA;
    double velocityDec;
};


#endif

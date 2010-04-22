// -*- LSST-C++ -*-
/******************************************************
 *
 * Author: Matthew Cleveland
 * File: findTracklets.h, cc
 * Purpose: To find tracklets.
 *
 ******************************************************/

#ifndef __FINDTRACKLETS_H__
#define __FINDTRACKLETS_H__

#include <vector>

#include "lsst/mops/Tracklet.h"
#include "lsst/mops/Detection.h"


namespace lsst {
    namespace mops {

/*****************************************************************
 * Main function
 *****************************************************************/
std::vector <Tracklet> 
findTracklets(const std::vector<Detection> &allDetections, 
	      double maxVelocity, double minVelocity);

    }} // close lsst::mops

#endif

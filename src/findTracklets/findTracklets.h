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

#include "../Tracklet.h"
#include "../Detection.h"


/*****************************************************************
 * Main function
 *****************************************************************/
std::vector <Tracklet> 
findTracklets(const std::vector<Detection> &allDetections, 
	      double maxVelocity, double minVelocity);

#endif

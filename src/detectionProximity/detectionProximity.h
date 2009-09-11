// -*- LSST-C++ -*-
/*
 *
 *
 */
#ifndef _DETECTION_PROXIMITY_
#define _DETECTION_PROXIMITY_


#include <stdlib.h>
#include <iostream>
#include <string>
#include <getopt.h>
#include <vector>
#include <fstream>
#include <utility> //for 'pair'

#include "../KDTree.h"
#include "../Detection.h"
#include "../collapseTrackletsAndPostfilters/TrackletCollapser.h"


/*
 * for each query point, find data points sufficiently nearby.
 * 
 * returns a vector of pairs of similar points; each pair has as its
 * first part an index into queryPoints and as its second part an
 * index into dataPoints.
 */
std::vector<std::pair <unsigned int, unsigned int> > 
detectionProximity(const std::vector<Detection>& queryPoints,
		   const std::vector<Detection>& dataPoints,
		   double distanceThreshold,
		   double brightnessThreshold,
		   double timeThreshold);

#endif

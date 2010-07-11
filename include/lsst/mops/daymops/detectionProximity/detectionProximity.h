// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
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

#include "lsst/mops/KDTree.h"
#include "lsst/mops/MopsDetection.h"


namespace lsst {
    namespace mops {

/*
 * for each query point, find data points sufficiently nearby.
 * 
 * returns a vector of pairs of similar points; each pair has as its
 * first part an index into queryPoints and as its second part an
 * index into dataPoints.
 */
std::vector<std::pair <unsigned int, unsigned int> > 
detectionProximity(const std::vector<MopsDetection>& queryPoints,
		   const std::vector<MopsDetection>& dataPoints,
                   double distanceThreshold,
		   double timeThreshold);

    }} // close lsst::mops

#endif

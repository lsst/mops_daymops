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
#ifndef _FIELD_PROXIMITY_
#define _FIELD_PROXIMITY_


#include <stdlib.h>
#include <iostream>
#include <string>
#include <getopt.h>
#include <vector>
#include <fstream>


#include "lsst/mops/KDTree.h"
#include "lsst/mops/daymops/fieldProximity/Field.h"
#include "lsst/mops/daymops/fieldProximity/TrackForFieldProximity.h"

namespace lsst {
    namespace mops {

/* returns a vector of pairs. for each pair: 

   pair.first should be an index into
   allTracks.  pair.second should be an index into queryFields.
   allTracks[pair.first] intersects queryFields[pair.second] 
*/

std::vector<std::pair <unsigned int, unsigned int>  >
fieldProximity(std::vector<FieldProximityTrack> allTracks,
               std::vector<Field> queryFields,
               double distThresh);


    }} // close lsst::mops

#endif

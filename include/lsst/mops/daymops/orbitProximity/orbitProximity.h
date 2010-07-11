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
#ifndef _ORBIT_PROXIMITY_
#define _ORBIT_PROXIMITY_


#include <iostream>
#include <string>
#include <getopt.h>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <utility>



#include "lsst/mops/Orbit.h"

namespace lsst {
    namespace mops {

/* orbitProximity: find similar orbits, where "similar" is
 * defined by the *Tolerance arguments. 
 * 
 * returns a vector of pairs. for each pair in pairs:
 * 
 *   dataOrbits[pair.first] is similar to 
 *   queryOrbits[pair.second]
 *
 */
std::vector<std::pair<unsigned int, unsigned int> > 
orbitProximity(std::vector<Orbit> dataOrbits, 
               std::vector<Orbit> queryOrbits,
               double perihelionTolerance,
               double eccentricityTolerance,
               double inclinationTolerance,
               double perihelionArgTolerance,
               double longitudeArgTolerance,
               double perihelionTimeTolerance);

    }} // close lsst::mops

#endif

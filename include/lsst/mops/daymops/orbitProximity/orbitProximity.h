// -*- LSST-C++ -*-
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

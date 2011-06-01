// -*- LSST-C++ -*-
/*
 *
 *
 */
#ifndef _FIELD_PROXIMITY_
#define _FIELD_PROXIMITY_


#include <stdlib.h>
#include <utility>
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

// queryFields may be modified - we will sort it by obs time.
void fieldProximity(const std::vector<FieldProximityTrack> &allTracks,
                    std::vector<Field> &queryFields,
                    std::vector<std::pair<unsigned int, unsigned int>  > &results,
                    double distThresh);

// legacy interface - will be slower because we copy the output
// vector, but needed to get unit tests compiling
std::vector<std::pair<unsigned int, unsigned int>  > 
fieldProximity(const std::vector<FieldProximityTrack> &allTracks,
               std::vector<Field> &queryFields,
               double distThresh);


    }} // close lsst::mops

#endif

// -*- LSST-C++ -*-
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


#include "../KDTree.h"
#include "../collapseTrackletsAndPostfilters/TrackletCollapser.h"
#include "Field.h"
#include "TrackForFieldProximity.h"

/* returns a vector of pairs. for each pair: 

   pair.first should be an index into
   allTracks.  pair.second should be an index into queryFields.
   allTracks[pair.first] intersects queryFields[pair.second] 
*/

std::vector<std::pair <unsigned int, unsigned int>  >
fieldProximity(std::vector<FieldProximityTrack> allTracks,
               std::vector<Field> queryFields,
               double distThresh);



#endif

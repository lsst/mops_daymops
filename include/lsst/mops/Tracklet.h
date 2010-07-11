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

namespace lsst {
namespace mops {

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


} } //close namespace lsst::mops

#endif

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
 * jmyers 7/29/08
 * 
 * collapse-tracklets is my specialized Hough transform code.  It is used for
 * quickly finding tracklets which describe similar linear motion and joining
 * them into longer tracklets.  Currently it is intended to be run from the
 * command-line, but the interface exported here should be fairly easy to port
 * into an LSST pipeline.
 *
 */



#ifndef LSST_COLLAPSE_TRACKLETS_H
#define LSST_COLLAPSE_TRACKLETS_H

#include <set>
#include <string>
#include <vector>

#include "lsst/mops/KDTree.h"
#include "lsst/mops/MopsDetection.h" 
#include "lsst/mops/Exceptions.h"
#include "lsst/mops/Tracklet.h"

namespace lsst {
    namespace mops {


    class TrackletCollapser {
    public:
                
        
        /* add all det indices from t1 (which are not already in t2) into t2. 
         * set t1 and t2.isCollapsed() = True.
         */
        void collapse(Tracklet &t1, Tracklet &t2);
        
        
        /* *pairs is modified - the Tracklets will have isCollapsed set. collapsedPairs will
         * actual output data.  I.e. if pairs contains similar tracklets [1,2]  and [2,3]  they will be
         * marked as collapsed, and [1,2,3] will be added to the collapsedPairs vector.*/
        void doCollapsingPopulateOutputVector(
            const std::vector<MopsDetection> * detections, 
            std::vector<Tracklet> &pairs,
            std::vector<double> tolerances, 
            std::vector<Tracklet> &collapsedPairs,
            bool useMinimumRMS, bool useBestFit, 
            bool useRMSFilt, double maxRMS, bool beVerbose);
      
        void setPhysicalParamsVector(const std::vector<MopsDetection> *trackletDets,
                                     std::vector<double> &physicalParams,
                                     double normalTime);
            

        void populateTrackletsForTreeVector(const std::vector<MopsDetection> *detections,
                                            const std::vector<Tracklet> * tracklets,
                                            std::vector<PointAndValue <unsigned int> >
                                            &trackletsForTree);

            
    };
    
    }} // close lsst::mops

#endif

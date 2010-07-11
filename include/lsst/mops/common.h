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
 * this is a place for common math functions and handy declarations
 * used by the KDTree class. */

#ifndef LSST_KDTREE_COMMON_H
#define LSST_KDTREE_COMMON_H



#include <vector>
#include <cmath>

#include <string>

namespace lsst {
namespace mops {
        
        enum GeometryType {
            EUCLIDEAN,
            CIRCULAR_DEGREES,
            CIRCULAR_RADIANS,
            RA_DEGREES,
            DEC_DEGREES,
            RA_RADIANS,
            DEC_RADIANS
        };

        class Constants {
        public:
            double epsilon() { return 1.e-10; }; 
            double rad_to_deg()  { return 180./M_PI; };
            double deg_to_rad()  { return M_PI / 180.; };

        };

        /* O(n)-time implementation returns median value of fv */
        double fastMedian(std::vector<double> fv);

        void printDoubleVec(std::vector<double> fv);

        /* returns whether a, b are within epsilon of each other */
        bool areEqual(double a, double b);
        
        /* returns distance from p1 to p2, ASSUMES that p1 and p2 have
           size==dimensions*/
        double euclideanDistance(std::vector<double> p1, std::vector<double> p2, 
                                unsigned int dimensions);

        double distance1D(double a, double b, GeometryType geo);

        /* return the length of the shortest path from a to b, where a, b 
           are angles.*/
        double circularShortestPathLen_Deg(double a, double b);
        
        double circularShortestPathLen_Rad(double a, double b);

        /* units in DEGREES.  use this if you have no prior knowledge about your
         * two regions of angular space. Any region which is 180 degrees will be assumed to 
         * encompass the whole circle and returns TRUE no matter what.*/
        bool angularRegionsOverlapSafe(double a1, double a2, double b1, double b2);

       /*
         * returns true iff the line segments aLo aHi and b1b2 overlap on the 1D
         * space of type GeometryType. 
         * 
         * for a circular space, we ASSUME that aLo is less than aHi, and aLo,
         * aHi may span up to 360 degrees (i.e., we always assume that aLo, aHi
         * does *not* cross 0).
         *
         * b1b2 is understood to be the shortest path * between b1 and b2.
         * E.g.:
         *
         * if aLo = 10, aHi = 310, then the length of aLo-aHi is 300.
         * if b1 = 10, b2 = 310, then the length of b1 to b2 is 60.
         *
         * This is useful for KDTreeNode hyperRectangle search; we have the data
         * partitioned by its numerical size as though euclidean; the searches
         * must be searching regions less than 180 degrees.
         */
        bool regionsOverlap1D(double aLo, double aHi, double b1, double b2, 
                              GeometryType type);


        // euclidean only, and aLo MUST be < aHi, bLo MUST be < bHi
        bool regionsOverlap1D_unsafe(double aLo, double aHi, double bLo, double bHi) ;

        /* for deg is any double value, return a value along [0, 369) */
        double convertToStandardDegrees(double deg);


        /* returns the great-circle distance between (RA0, Dec0) and (RA1, Dec1) - 
         * that is, the shortest path along the surface of a unit sphere between
         * these two points (on the unit sphere). All units in degrees. */
        double angularDistanceRADec_deg(double RA0, double Dec0, double RA1, double Dec1);


        /* 
         * arcToRA: given a lattitude in declination (degrees) and an arc length
         * in degrees, return the number of degrees in RA which corresponding to
         * that arc length at the given declination.

         * E.G, if you are at RA,Dec = (100,0) and travel 1 degree to (101,0), you have
         * moved 1 degree in RA and 1 degree in 'actual' distance Ergo
         * arcToRA(0,1) == 1.  

         * however, very close to the north pole at RA,Dec = (0,89), you can walk halfway
         * 'around the earth' to (180,89) and only cross 2 degrees of arc. Ergo
         * arcToRA(89,2) == 180.
         */
        double arcToRA(double dec, double arcDeg);       
        
        double minOfTwo(double a, double b);

        double maxOfTwo(double a, double b);

        bool guessBoolFromStringOrGiveErr(std::string guessStr, std::string errStr);

        std::string boolToString(bool b);

        std::string stringToLower(std::string strToConvert);
        

}} // close namespace lsst::mops

#endif

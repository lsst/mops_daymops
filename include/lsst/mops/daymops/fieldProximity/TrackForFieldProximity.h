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
 */

#ifndef LSST_MOPS_FP_TRACK_H
#define LSST_MOPS_FP_TRACK_H

#include <string>


class FieldProximityPoint {
public: 
    void setRA(double w) {RA = w;}
    void setDec(double w) {dec = w;}
    void setEpochMJD(double w) {MJD = w;}

    double getEpochMJD() const {return MJD;}
    double getRA() const {return RA; }
    double getDec() const {return dec; }

private: 
    double MJD;

    //all in degrees
    double RA;
    double dec;
    double radius;
    
};

/* fieldProximity deals with a very simple notion of 'track' as a set of unique
 * points.  This is quite different from the tracks/tracklets we use elsewhere,
 * which are linkages between MopsDetections (DIASources).  Hence, fieldProximity
 * gets its own special track class
 */
class FieldProximityTrack {
public:

    FieldProximityTrack() { myID = ""; myPoints.clear(); }

    void setID(std::string s) { myID = s; }

    void setPoints(std::vector<FieldProximityPoint> w) { myPoints = w; }

    void addPoint(FieldProximityPoint w) { myPoints.push_back(w); }

    std::string getID() const { return myID; }
    std::vector<FieldProximityPoint> getPoints() const { return myPoints; }

private:
    std::string myID;
    std::vector<FieldProximityPoint> myPoints;
};



#endif

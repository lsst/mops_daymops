// -*- LSST-C++ -*-

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

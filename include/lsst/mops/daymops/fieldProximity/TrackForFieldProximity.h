// -*- LSST-C++ -*-

/*
 * jmyers 7/29/08
 * 
 */

#ifndef LSST_MOPS_FP_TRACK_H
#define LSST_MOPS_FP_TRACK_H

#include <string>


namespace lsst { namespace mops {

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

    FieldProximityTrack() { myID = 0; myPoints.clear(); }

    void setID(uint s) { myID = s; }

    void setPoints(std::vector<FieldProximityPoint> w) { myPoints = w; }

    void addPoint(FieldProximityPoint w) { myPoints.push_back(w); }

    unsigned int getID() const { return myID; }
    const std::vector<FieldProximityPoint> * getPoints() const { return &myPoints; }

private:
    unsigned int myID;
    std::vector<FieldProximityPoint> myPoints;
};

}} // close lsst::mops

#endif

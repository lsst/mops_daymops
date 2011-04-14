// -*- LSST-C++ -*-

/*
 * jmyers 7/29/08
 *
 * "lightened up" on Apr 23 2010
 *
 * This is the "lightened" Detection class; it contains only the data needed by
 * DayMops to do its work.  It does not attempt to hold as much data as real
 * DiaSource, or even as would be held in Miti format.  But since we expect data
 * transmission to be an issue in LSST, we're cutting down on the size of data
 * here.
 * 
 */

#ifndef LSST_MOPS_DETECTION_H
#define LSST_MOPS_DETECTION_H


#include <fstream>
#include <iostream>
#include <string>


namespace lsst {
namespace mops {

class MopsDetection {
public:

    //create an empty detection
    MopsDetection();

    //create a "lightweight" detection.
    MopsDetection(long int ID, double epochMJD, double RA, double Dec, double RaErr=0, double DecErr=0);
    
    // the MITI format holds lots more data than ID, MJD, Ra, and Dec.  
    // We don't store those items, though!
    void fromMITIString(std::string);

    //getters...  note that if you don't set values, you'll get bogus
    // placeholder values from the "getters" - they used to send exceptions in
    // this case, but to save the extra byte, I dropped the "initialized" flag.
    long int getID() const ;
    double getEpochMJD() const ;
    double getRA() const ;
    double getDec() const ;
    double getRaErr() const ;
    double getDecErr() const ;
    double getRaTopoCorr() const ;

    void setID(long int newId);
    void setEpochMJD(double newMjd);
    void setRA(double newRa);
    void setDec(double newDec);
    void setRaErr(double RaErr);
    void setDecErr(double DecErr);
    void calculateTopoCorr(double obsLat, double obsLong);

private:

    long int ID;
    double MJD;
    double RA;
    double dec;
    double RaErr;
    double DecErr;
    double RaTopoCorr;
};

}} // close namespace lsst::mops

#endif

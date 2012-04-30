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
    MopsDetection(long int ID, double epochMJD, double RA, double Dec, double RaErr=0, double DecErr=0, int ssmId=-1, long int obsHistId=-1, double snr=-1, double mag=-1);
    
    void fromString(std::string);

    //getters...  note that if you don't set values, you'll get bogus
    // placeholder values from the "getters" - they used to send exceptions in
    // this case, but to save the extra byte, I dropped the "initialized" flag.
    long int getID() const ;
    long int getImageID() const;
    double getEpochMJD() const ;
    double getRA() const ;
    double getDec() const ;
    double getMag() const ; 
    double getSNR() const ; 
    int getSsmId() const;
    double getRaErr() const ;
    double getDecErr() const ;
    double getRaTopoCorr() const ;

    void setID(long int newId);
    void setImageID(long int);
    void setEpochMJD(double newMjd);
    void setRA(double newRa);
    void setDec(double newDec);
    void setMag(double) ; 
    void setSNR(double) ; 
    void setSsmId(int newSsmId);
    void setRaErr(double RaErr);
    void setDecErr(double DecErr);
    void calculateTopoCorr();

    static void setObservatoryLocation(double obsLat, double obsLong);

private:

    static double obsLat;
    static double obsLong;

    long int ID;
    double MJD;
    double RA;
    double dec;
    double RaErr;
    double DecErr;
    double RaTopoCorr;
    int ssmId;
    long int imageID;
    double mag;
    double snr;
};

}} // close namespace lsst::mops

#endif

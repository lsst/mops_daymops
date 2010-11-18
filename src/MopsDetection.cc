// -*- LSST-C++ -*-
#include <iomanip>
#include <sstream>

#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/Exceptions.h"



/*
 * jmyers 7/29/08
 * 
 */

namespace lsst {
    namespace mops {

MopsDetection::MopsDetection()
{
    RA = -380;
    dec = -380;
    MJD = 0;
}


MopsDetection::MopsDetection(long int ID, double epochMJD, double RA, double Dec) 
{
    this->ID = ID;
    MJD = epochMJD;
    this->RA = RA;
    this->dec = Dec;
}



void MopsDetection::setID(long int newId)
{
    ID = newId;
}


void MopsDetection::setEpochMJD(double newMjd)
{
    MJD = newMjd;
}

void MopsDetection::setRA(double newRa)
{
    RA = newRa;
}


void MopsDetection::setDec(double newDec)
{
    dec = newDec;
}
        

long int MopsDetection::getID() const 
{
    return ID;

}

double MopsDetection::getEpochMJD() const  
{
    return MJD;
}


double MopsDetection::getRA()  const 
{
    return RA;
}


double MopsDetection::getDec()  const 
{
    return dec;
}




void MopsDetection::fromMITIString(std::string mitiString) {
    /*     --------- PANSTARRS Input File Format (from Larry Denneau's Spec):
           ID EPOCH_MJD RA_DEG DEC_DEG MAG OBSCODE OBJECT_NAME LENGTH ANGLE [ETIME]*/
    // read these parameters to local vars; we don't save them.
    double mag;
    std::string obscode;
    std::string objName;
    double length;
    double angle;
    double etime;
    bool hasETime;
    std::istringstream ss(mitiString);
    /* make SS raise an exception if reading doesn't happen correctly */
    ss.exceptions(std::ifstream::failbit | std::ifstream::badbit);    
    try {
        ss >> ID;
        ss >> MJD;
        ss >> RA;
        ss >> dec;
        ss >> mag;
        ss >> obscode;
        ss >> objName;
        ss >> length;
        ss >> angle;
    }
    catch (...) {
        throw LSST_EXCEPT(BadParameterException, 
                          "Badly-formatted MITI string\n");
    }
    try {
        ss >> etime;
        hasETime = true;
    }
    catch (...) {
        hasETime = false;
    }
}





} } // close lsst::mops namespace

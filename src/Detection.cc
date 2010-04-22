// -*- LSST-C++ -*-
#include <iomanip>
#include <sstream>

#include "lsst/mops/Detection.h"
#include "lsst/mops/Exceptions.h"



/*
 * jmyers 7/29/08
 * 
 */

namespace lsst {
    namespace mops {

Detection::Detection()
{
    hasETime = false;
    initialized = false;
    fileIndex = -1;
}


Detection::Detection(long int ID, double epochMJD, double RA, double Dec, 
                     int obsCode, std::string objName, double mag, 
                     double elongationLength, double elongationAngle)
{
    hasETime = false;
    initialized = true;
    this->ID = ID; 
    MJD = epochMJD;
    this->RA = RA;
    this->dec = Dec;
    this->obscode = obsCode;
    this->objName = objName;
    this->mag = mag;
    length = elongationLength;
    angle = elongationAngle;
    fileIndex = -1;
}

Detection::Detection(long int ID, double epochMJD, double RA, double Dec, 
                     int obsCode, std::string objName, double mag, 
                     double elongationLength, double elongationAngle,
                     double exposureTime)
{
    hasETime = true;
    initialized = true;
    this->ID = ID; 
    MJD = epochMJD;
    this->RA = RA;
    this->dec = Dec;
    this->obscode = obsCode;
    this->objName = objName;
    this->mag = mag;
    length = elongationLength;
    angle = elongationAngle;
    etime = exposureTime;
    fileIndex = -1;
}



Detection::Detection(long int ID, double epochMJD, double RA, double Dec) 
{
    hasETime = false;
    initialized = true;
    this->ID = ID;
    MJD = epochMJD;
    this->RA = RA;
    this->dec = Dec;
    this->obscode = 0;
    this->mag = 0.;
    fileIndex = -1;
}



long int Detection::getID() const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "Detection: request for ID from uninitialized Detection\n");
        return -1;
    }
    else {
        return ID;
    }

}

double Detection::getEpochMJD() const  
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "Detection: request for EpochMJD from uninitialized Detection\n");
        return -1;
    }
    else {
        return MJD;
    }
}


double Detection::getRA()  const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "Detection: request for RA from uninitialized Detection\n");
        return -1;
    }
    else {
        return RA;
    }
}


double Detection::getDec()  const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "Detection: request for Dec from uninitialized Detection\n");
        return -1;
    }
    else {
        return dec;
    }
}


int  Detection::getObscode() const  
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "Detection: request for ObsCode from uninitialized Detection\n");
        return -1;
    }
    else {
        return obscode;
    }

}

std::string Detection::getObjName() const  
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "Detection: request for ObjName from uninitialized Detection\n");
        return "";
    }
    else {
        return objName;
    }
}

double Detection::getMag()  const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "Detection: request for magnitude from uninitialized Detection\n");
        return -1000;
    }
    else {
        return mag;
    }
}


double Detection::getLength()  const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "Detection: request for length from uninitialized Detection\n");
        return -1000;
    }
    else {
        return length;
    }
}



double Detection::getAngle() const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "Detection: request for angle from uninitialized Detection\n");
        return -1000;
    }
    else {
        return angle;
    }
    
}


bool Detection::hasExposureTime() const {
    return hasETime;
}

double Detection::getExposureTime() const 
{
    if ((!initialized) || (!hasETime)) {
        throw LSST_EXCEPT(UninitializedException, 
                          "Detection: request for exposure from Detection which is uninitialized or has no etime data\n");
        return -1000;
    }
    else {
        return etime;
    }
    
}




void Detection::fromMITIString(std::string mitiString) {
    /*     --------- PANSTARRS Input File Format (from Larry Denneau's Spec):
           ID EPOCH_MJD RA_DEG DEC_DEG MAG OBSCODE OBJECT_NAME LENGTH ANGLE [ETIME]*/
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
    initialized = true;
}



void Detection::setFileIndex(int v)
{
  fileIndex = v;
}



int Detection::getFileIndex() const
{

  if ((!initialized)) {
    throw LSST_EXCEPT(UninitializedException, 
		      "Detection: request for exposure from Detection which is uninitialized or has no etime data\n");
    return -1000;
  }
  else {
    return fileIndex;
  }
}



    } } // close lsst::mops namespace

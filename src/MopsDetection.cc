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
    hasETime = false;
    initialized = false;
    fileIndex = -1;
}


MopsDetection::MopsDetection(long int ID, double epochMJD, double RA, double Dec, 
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

MopsDetection::MopsDetection(long int ID, double epochMJD, double RA, double Dec, 
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



MopsDetection::MopsDetection(long int ID, double epochMJD, double RA, double Dec) 
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



long int MopsDetection::getID() const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MopsDetection: request for ID from uninitialized MopsDetection\n");
        return -1;
    }
    else {
        return ID;
    }

}

double MopsDetection::getEpochMJD() const  
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MopsDetection: request for EpochMJD from uninitialized MopsDetection\n");
        return -1;
    }
    else {
        return MJD;
    }
}


double MopsDetection::getRA()  const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MopsDetection: request for RA from uninitialized MopsDetection\n");
        return -1;
    }
    else {
        return RA;
    }
}


double MopsDetection::getDec()  const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MopsDetection: request for Dec from uninitialized MopsDetection\n");
        return -1;
    }
    else {
        return dec;
    }
}


int  MopsDetection::getObscode() const  
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MopsDetection: request for ObsCode from uninitialized MopsDetection\n");
        return -1;
    }
    else {
        return obscode;
    }

}

std::string MopsDetection::getObjName() const  
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MopsDetection: request for ObjName from uninitialized MopsDetection\n");
        return "";
    }
    else {
        return objName;
    }
}

double MopsDetection::getMag()  const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MopsDetection: request for magnitude from uninitialized MopsDetection\n");
        return -1000;
    }
    else {
        return mag;
    }
}


double MopsDetection::getLength()  const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MopsDetection: request for length from uninitialized MopsDetection\n");
        return -1000;
    }
    else {
        return length;
    }
}



double MopsDetection::getAngle() const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MopsDetection: request for angle from uninitialized MopsDetection\n");
        return -1000;
    }
    else {
        return angle;
    }
    
}


bool MopsDetection::hasExposureTime() const {
    return hasETime;
}

double MopsDetection::getExposureTime() const 
{
    if ((!initialized) || (!hasETime)) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MopsDetection: request for exposure from MopsDetection which is uninitialized or has no etime data\n");
        return -1000;
    }
    else {
        return etime;
    }
    
}




void MopsDetection::fromMITIString(std::string mitiString) {
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



void MopsDetection::setFileIndex(int v)
{
  fileIndex = v;
}



int MopsDetection::getFileIndex() const
{

  if ((!initialized)) {
    throw LSST_EXCEPT(UninitializedException, 
		      "MopsDetection: request for exposure from MopsDetection which is uninitialized or has no etime data\n");
    return -1000;
  }
  else {
    return fileIndex;
  }
}



    } } // close lsst::mops namespace

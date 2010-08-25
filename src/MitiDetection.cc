// -*- LSST-C++ -*-
#include <iomanip>
#include <sstream>

#include "lsst/mops/MitiDetection.h"
#include "lsst/mops/Exceptions.h"



/*
 * jmyers 7/29/08
 * 
 */

namespace lsst {
    namespace mops {

MitiDetection::MitiDetection()
{
    hasETime = false;
    initialized = false;
    fileIndex = -1;
}


MitiDetection::MitiDetection(long int ID, double epochMJD, double RA, double Dec, 
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

MitiDetection::MitiDetection(long int ID, double epochMJD, double RA, double Dec, 
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






long int MitiDetection::getID() const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MitiDetection: request for ID from uninitialized MitiDetection\n");
        return -1;
    }
    else {
        return ID;
    }

}

double MitiDetection::getEpochMJD() const  
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MitiDetection: request for EpochMJD from uninitialized MitiDetection\n");
        return -1;
    }
    else {
        return MJD;
    }
}


double MitiDetection::getRA()  const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MitiDetection: request for RA from uninitialized MitiDetection\n");
        return -1;
    }
    else {
        return RA;
    }
}


double MitiDetection::getDec()  const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MitiDetection: request for Dec from uninitialized MitiDetection\n");
        return -1;
    }
    else {
        return dec;
    }
}


int  MitiDetection::getObscode() const  
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MitiDetection: request for ObsCode from uninitialized MitiDetection\n");
        return -1;
    }
    else {
        return obscode;
    }

}

std::string MitiDetection::getObjName() const  
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MitiDetection: request for ObjName from uninitialized MitiDetection\n");
        return "";
    }
    else {
        return objName;
    }
}

double MitiDetection::getMag()  const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MitiDetection: request for magnitude from uninitialized MitiDetection\n");
        return -1000;
    }
    else {
        return mag;
    }
}


double MitiDetection::getLength()  const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MitiDetection: request for length from uninitialized MitiDetection\n");
        return -1000;
    }
    else {
        return length;
    }
}



double MitiDetection::getAngle() const 
{
    if (!initialized) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MitiDetection: request for angle from uninitialized MitiDetection\n");
        return -1000;
    }
    else {
        return angle;
    }
    
}


bool MitiDetection::hasExposureTime() const {
    return hasETime;
}

double MitiDetection::getExposureTime() const 
{
    if ((!initialized) || (!hasETime)) {
        throw LSST_EXCEPT(UninitializedException, 
                          "MitiDetection: request for exposure from MitiDetection which is uninitialized or has no etime data\n");
        return -1000;
    }
    else {
        return etime;
    }
    
}




void MitiDetection::fromMITIString(std::string mitiString) {
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



void MitiDetection::setFileIndex(int v)
{
  fileIndex = v;
}



int MitiDetection::getFileIndex() const
{

  if ((!initialized)) {
    throw LSST_EXCEPT(UninitializedException, 
		      "MitiDetection: request for exposure from MitiDetection which is uninitialized or has no etime data\n");
    return -1000;
  }
  else {
    return fileIndex;
  }
}



    } } // close lsst::mops namespace

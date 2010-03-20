// -*- LSST-C++ -*-

/*
 * jmyers 7/29/08
 * 
 */

#ifndef LSST_MOPS_DETECTION_H
#define LSST_MOPS_DETECTION_H


#include <fstream>
#include <iostream>
#include <string>



class Detection {
public:

    //create an empty detection
    Detection();


    // these two added 6/19/09 (jmyers)
    // create a detection with initialized data, but no exposure time
    Detection(long int ID, double epochMJD, double RA, double Dec, 
              std::string obsCode, std::string objName, double mag,
              double elongationLength, double elongationAngle);

    //create a detection with initialized data, including exposure time data
    Detection(long int ID, double epochMJD, double RA, double Dec, 
              std::string obsCode, std::string objName, double mag, 
              double elongationLength, double elongationAngle, 
              double exposureTime);

    // a detection is only initialized after we have assigned
    // it values from a string.  If you try to use a getter
    // function before data is initialized, you will get an
    // exception.
    bool isInitialized() const {return initialized; }    

    void fromMITIString(std::string);

    //getters...
    long int getID() const ;
    double getEpochMJD() const ;
    double getRA() const ;
    double getDec() const ;
    std::string  getObscode() const ;
    std::string getObjName() const ;
    double getMag() const ;
    double getLength() const ;
    double getAngle() const ;

    // exposure time is optional in PanSTARRS MITI format.
    // also optional in our constructors.
    bool hasExposureTime() const ;
    double getExposureTime() const ;
    void setFileIndex(int);//{fileIndex = v;};
    
    int getFileIndex() const;//{return fileIndex;};

private:

    int fileIndex;
    long int ID;
    double MJD;
    double RA;
    double dec;
    std::string obscode;
    double mag;
    std::string objName;
    double length;
    double angle;
    double etime;
    bool hasETime;
    bool initialized;
};



#endif

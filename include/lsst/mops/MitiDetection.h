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
 * jmyers 4/23/10
 *
 * This was formerly the "overweight" version of what is now called MopsDetection.
 * 
 * It uses more information than is strictly needed for the DayMops apps - but
 * it corresponds nicely to PanSTARRS MITI format.  It inherits from MopsDetection
 * and adds lots of new fields.
 * 
 * If memory is an issue - and it probably is - use MopsDetection.  If it's not
 * an issue, consider LSST DiaSource.  If you're not LSST and don't want to use
 * LSST software, use this.
 */

#ifndef LSST_MITI_DETECTION_H
#define LSST_MITI_DETECTION_H


#include <fstream>
#include <iostream>
#include <string>

#include "lsst/mops/MopsDetection.h"

namespace lsst {
namespace mops {

    class MitiDetection : public MopsDetection {
public:

    //create an empty detection
    MitiDetection();


    // these two added 6/19/09 (jmyers)
    // create a detection with initialized data, but no exposure time
    MitiDetection(long int ID, double epochMJD, double RA, double Dec, 
              int obsCode, std::string objName, double mag,
              double elongationLength, double elongationAngle);

    //create a detection with initialized data, including exposure time data
    MitiDetection(long int ID, double epochMJD, double RA, double Dec, 
              int obsCode, std::string objName, double mag, 
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
    int  getObscode() const ;
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
    int obscode;
    double mag;
    std::string objName;
    double length;
    double angle;
    double etime;
    bool hasETime;
    bool initialized;
};

}} // close namespace lsst::mops

#endif

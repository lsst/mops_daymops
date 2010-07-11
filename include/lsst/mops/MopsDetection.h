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
    MopsDetection(long int ID, double epochMJD, double RA, double Dec);
    
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

private:

    long int ID;
    double MJD;
    double RA;
    double dec;
};

}} // close namespace lsst::mops

#endif

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
 */

#ifndef LSST_MOPS_FIELD_H
#define LSST_MOPS_FIELD_H

#include <string>

class Field {
public:

    void setFileIndex(int);
    void setFieldID(std::string);
    void setEpochMJD(double);
    void setRA(double);
    void setDec(double);
    void setRadius(double);

    int getFileIndex() const;
    std::string getFieldID() const;
    double getEpochMJD() const;
    double getRA() const;
    double getDec() const;
    double getRadius() const;

private:

    int fileIndex;
    std::string fieldID;

    double MJD;

    //all in degrees
    double RA;
    double dec;
    double radius;
};



#endif

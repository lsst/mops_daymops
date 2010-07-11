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
 *
 */

#include "lsst/mops/daymops/fieldProximity/Field.h"
#include "lsst/mops/Exceptions.h"


/**************************************************
 * SETTERS
 **************************************************/
void Field::setFileIndex(int f)
{
    fileIndex = f;
}

void Field::setFieldID(std::string f)
{
    fieldID = f;
}

void Field::setEpochMJD(double e)
{
    MJD = e;
}

void Field::setRA(double ra)
{
    RA = ra;
}

void Field::setDec(double d)
{
    dec = d;
}

void Field::setRadius(double r)
{
    radius = r;
}



/**************************************************
 * GETTERS
 **************************************************/
int Field::getFileIndex() const
{
    return fileIndex;
}

std::string Field::getFieldID() const
{
    return fieldID;
}

double Field::getEpochMJD() const
{
    return MJD;
}

double Field::getRA() const
{
    return RA;
}

double Field::getDec() const
{
    return dec;
}

double Field::getRadius() const
{
    return radius;
}

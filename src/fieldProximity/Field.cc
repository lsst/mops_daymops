// -*- LSST-C++ -*-
/*
 *
 */

#include "Field.h"
#include "../Exceptions.h"


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

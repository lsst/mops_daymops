// -*- LSST-C++ -*-
/*
 *
 */

#include "lsst/mops/daymops/fieldProximity/Field.h"
#include "lsst/mops/Exceptions.h"

namespace lsst { namespace mops {

/**************************************************
 * SETTERS
 **************************************************/

void Field::setFieldID(unsigned int f)
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

unsigned int Field::getFieldID() const
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


}} // close lsst::mops

// -*- LSST-C++ -*-

/*
 * jmyers 7/29/08
 * 
 */

#ifndef LSST_MOPS_FIELD_H
#define LSST_MOPS_FIELD_H

#include <string>

namespace lsst { namespace mops {

class Field {
public:

    void setFieldID(unsigned int);
    void setEpochMJD(double);
    void setRA(double);
    void setDec(double);
    void setRadius(double);

    unsigned int getFieldID() const;
    double getEpochMJD() const;
    double getRA() const;
    double getDec() const;
    double getRadius() const;

private:

    unsigned int fieldID;

    double MJD;

    //all in degrees
    double RA;
    double dec;
    double radius;
};

}} // close lsst::mops

#endif

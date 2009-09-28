// -*- LSST-C++ -*-

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

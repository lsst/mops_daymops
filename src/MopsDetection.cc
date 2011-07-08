// -*- LSST-C++ -*-
#include <iomanip>
#include <sstream>

#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/Exceptions.h"

#include "slalib.h"
#include "slamac.h"

#undef DEBUG

/*
 * jmyers 7/29/08
 * 
 */

namespace lsst {
    namespace mops {

double MopsDetection::obsLat;
double MopsDetection::obsLong;

MopsDetection::MopsDetection()
{
    RA = -380;
    dec = -380;
    MJD = 0;
}


 MopsDetection::MopsDetection(long int ID, double epochMJD, double RA, double Dec, double RaErr, double DecErr, int ssmId) 
{
    this->ID = ID;
    MJD = epochMJD;
    this->RA = RA;
    this->dec = Dec;
    this->RaErr = RaErr;
    this->DecErr = DecErr;
    this->ssmId = ssmId;
}



void MopsDetection::setID(long int newId)
{
    ID = newId;
}

void MopsDetection::setSsmId(int newSsmId)
{
    ssmId = newSsmId;
}


void MopsDetection::setEpochMJD(double newMjd)
{
    MJD = newMjd;
}

void MopsDetection::setRA(double newRa)
{
    RA = newRa;
}


void MopsDetection::setDec(double newDec)
{
    dec = newDec;
}
        
void MopsDetection::setRaErr(double newRaErr)
{
    RaErr = newRaErr;
}
        
void MopsDetection::setDecErr(double newDecErr)
{
    DecErr = newDecErr;
}
 
void MopsDetection::setObservatoryLocation(double lat, double longitude)
{
    obsLat = lat;
    obsLong = longitude;
}

long int MopsDetection::getID() const 
{
    return ID;

}

int MopsDetection::getSsmId() const 
{
    return ssmId;

}

double MopsDetection::getEpochMJD() const  
{
    return MJD;
}


double MopsDetection::getRA()  const 
{
    return RA;
}


double MopsDetection::getDec()  const 
{
    return dec;
}

double MopsDetection::getRaErr()  const 
{
    return RaErr;
}

double MopsDetection::getDecErr()  const 
{
    return DecErr;
}

double MopsDetection::getRaTopoCorr()  const 
{
    return RaTopoCorr ;
}




void MopsDetection::fromMITIString(std::string mitiString) {
    /*     --------- PANSTARRS Input File Format (from Larry Denneau's Spec):
           ID EPOCH_MJD RA_DEG DEC_DEG MAG OBSCODE OBJECT_NAME LENGTH ANGLE [ETIME]*/
    // read these parameters to local vars; we don't save them.
    double mag;
    std::string obscode;
    double length;
    double angle;
    double etime;
    bool hasETime;
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
        ss >> ssmId;
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
}


void MopsDetection::calculateTopoCorr() {

    double obsLatRad, obsLongRad;

    obsLatRad = obsLat*DD2R;
    obsLongRad = obsLong*DD2R;
    
    double localAppSidTime = slaGmst(MJD - slaDt(slaEpj(MJD))/86400.0) + obsLongRad;

    double geoPosVel[6]; // observing position (and velocity) in AU, AU/sec
    slaPvobs(obsLatRad, 0, localAppSidTime, geoPosVel);
    
    double raRad, decRad;
    raRad = RA*DD2R;
    decRad = dec*DD2R;
    
    float rho[3];   // geocentric unit 3-vector to object
    slaCs2c(raRad, decRad, rho);

    // add geoPos to rho (multiplied by 1 AU) to get the topocentric vector to the object
    float rhoTopo[3];
    rhoTopo[0] = rho[0] + geoPosVel[0];
    rhoTopo[1] = rho[1] + geoPosVel[1];
    rhoTopo[2] = rho[2] + geoPosVel[2];

    // calculate the topocentric ra, dec

    float raTopo, decTopo;
    slaCc2s(rhoTopo, &raTopo, &decTopo);

    double deltaRa = raTopo - raRad;

    // make sure result is in right quadrant

    if (deltaRa < -DPIBY2) {
        deltaRa += D2PI;
    } else if (deltaRa > DPIBY2) {
        deltaRa -= D2PI;
    }

    RaTopoCorr = deltaRa*DR2D;

#ifdef DEBUG
    std::cerr << "topo_corr: " << MJD << ' ' << RA << ' ' << localAppSidTime << ' ' << RaTopoCorr << '\n';
#endif
    
}


} } // close lsst::mops namespace

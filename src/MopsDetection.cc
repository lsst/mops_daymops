// -*- LSST-C++ -*-
#include <iomanip>
#include <sstream>

#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/Exceptions.h"

extern "C" {
#include "pal.h"
}
#include "palmac.h"

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



MopsDetection::MopsDetection(long int ID, double epochMJD, double RA, double Dec, double RaErr, double DecErr, int ssmId, long int obsHistId, double snr, double mag) 
{
    this->ID = ID;
    MJD = epochMJD;
    this->RA = RA;
    this->dec = Dec;
    this->RaErr = RaErr;
    this->DecErr = DecErr;
    this->ssmId = ssmId;
    this->imageID = obsHistId;
    this->snr = snr;
    this->mag = mag;
}



void MopsDetection::setID(long int newId)
{
    ID = newId;
}

void MopsDetection::setImageID(long int newiid) 
{
    imageID = newiid;
}

void MopsDetection::setMag(double newMag) 
{
    mag = newMag;
}

void MopsDetection::setSNR(double newSnr) 
{
    snr=newSnr;
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

long int MopsDetection::getImageID() const
{
    return imageID;
}

double MopsDetection::getMag() const
{
    return mag;
}

double MopsDetection::getSNR() const
{
    return snr;
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




void MopsDetection::fromString(std::string diaStr) {
    // fullerDiaSource format:

    //  diaId obsHistId/imageId ssmId RA Dec MJD mag SNR

    std::istringstream ss(diaStr);
    /* make SS raise an exception if reading doesn't happen correctly */
    ss.exceptions(std::ifstream::failbit | std::ifstream::badbit);    
    try {
        ss >> ID;
        ss >> imageID;
        ss >> ssmId;
        ss >> RA;
        ss >> dec;
        ss >> MJD;
        ss >> mag;
        ss >> snr;
    }
    catch (...) {
        throw LSST_EXCEPT(BadParameterException, 
                          "Badly-formatted DiaSource string. Note that MITI is no longer supported\n");
    }
}


void MopsDetection::calculateTopoCorr() {

    double obsLatRad, obsLongRad;

    obsLatRad = obsLat*DD2R;
    obsLongRad = obsLong*DD2R;
    
    double localAppSidTime = palGmst(MJD - slaDt(palEpj(MJD))/86400.0) + obsLongRad;

    double geoPosVel[6]; // observing position (and velocity) in AU, AU/sec
    palPvobs(obsLatRad, 0, localAppSidTime, geoPosVel);
    
    double raRad, decRad;
    raRad = RA*DD2R;
    decRad = dec*DD2R;
    
    float rho[3];   // geocentric unit 3-vector to object
    palCs2c(raRad, decRad, rho);

    // add geoPos to rho (multiplied by 1 AU) to get the topocentric vector to the object
    float rhoTopo[3];
    rhoTopo[0] = rho[0] + geoPosVel[0];
    rhoTopo[1] = rho[1] + geoPosVel[1];
    rhoTopo[2] = rho[2] + geoPosVel[2];

    // calculate the topocentric ra, dec

    float raTopo, decTopo;
    palCc2s(rhoTopo, &raTopo, &decTopo);

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

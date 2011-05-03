#include <algorithm>
#include <set>
#include <iostream>

#include "lsst/mops/Track.h"
#include "lsst/mops/Exceptions.h"

#undef DEBUG

namespace lsst { namespace mops {

Track::Track()
{
    epoch = 0;
}


void Track::addTracklet(unsigned int trackletIndex, 
			const Linkage *t, 
			const std::vector<MopsDetection> & allDets)
{
     componentTrackletIndices.insert(trackletIndex);
     std::set<unsigned int> tDets = t->getComponentDetectionIndices();
     std::cout << " inserting a new track(let) to this track. \n"
	       << " initial size is " << componentDetectionIndices.size()
	       << " \n";
     std::set<unsigned int>::const_iterator tDetIter;
	  for (tDetIter = tDets.begin();
	       tDetIter != tDets.end();
	       tDetIter++) {
	  componentDetectionIndices.insert(*tDetIter);
	  componentDetectionDiaIds.insert(allDets.at(*tDetIter).getID());
     }
     std::cout << " after insert size is "
	       << componentDetectionIndices.size()
	       << " \n";

}
	  




void Track::addDetection(unsigned int detIndex, const std::vector<MopsDetection> & allDets)
{
    // TDB: might be nice to check that we haven't been given two dets
    // from the same image time.  this will enforce some safety with
    // quadratic fitting. not critical for now, since linkTracklets
    // never does anything that dumb.
    componentDetectionIndices.insert(detIndex);
    componentDetectionDiaIds.insert(allDets.at(detIndex).getID());
}
	  
const std::set<unsigned int> Track::getComponentDetectionIndices() const
{
    const std::set<unsigned int>copy(componentDetectionIndices);
    return copy;
}

const std::set<unsigned int> Track::getComponentDetectionDiaIds() const
{
    const std::set<unsigned int>copy(componentDetectionDiaIds);
    return copy;
}


Track & Track::operator= (const Track &other) {
    componentTrackletIndices = other.componentTrackletIndices;
    componentDetectionIndices = other.componentDetectionIndices;
    componentDetectionDiaIds = other.componentDetectionDiaIds;
    raFunc = other.raFunc;
    decFunc = other.decFunc;
    chisqRa = other.chisqRa;
    chisqDec = other.chisqDec;
    probChisqRa = probChisqRa;
    probChisqDec = probChisqDec;
    epoch = other.epoch;
    myId = other.myId;
    numUniqueNights = other.numUniqueNights;
    return *this;
}





/* returns SSM ID of underlying object or -1 if a false track.*/
int Track::getObjectId(const std::vector<MopsDetection> &allDets) const
{
     std::set<unsigned int>::const_iterator detIter = componentDetectionIndices.begin();
     if (componentDetectionIndices.size() == 0) {
	  throw LSST_EXCEPT(UninitializedException,
			    "Object Id requested from empty track.");	  
     }
     // get a first guess at our ID. 
     int underlyingId = allDets.at(*detIter).getSsmId();
     // if any other ID is != underlyingId, then we're a false track.
     // note that if our initial underlying ID is -1 (noise) then
     // we'll return -1 no matter what.
     for (detIter = componentDetectionIndices.begin();
	  detIter != componentDetectionIndices.end();
	  detIter++) {
	  int thisId = allDets.at(*detIter).getSsmId();
	  if (thisId != underlyingId) {
	       return -1;
	  }
     }
     return underlyingId;

}



/* jmyers - these implement the interface inherited from Linkage; it
 * will be used by TrackletTree (which should probably be renamed) and
 * the accel comparison functions in linkTrackets.*/
void Track::calculateBestFitFunc(const std::vector<MopsDetection> &allDets)
{
     // is it wise to use the higher-order function for accel prune, etc...?
     // probably not - we're trying to use acceleration as a filter, 
     // higher-order factors could cause great confusion.
     calculateBestFitQuadratic(allDets, false);
}


void Track::getFitFunction(double& epoch, std::vector<double> &raFunc,
			   std::vector<double> &decFunc)
{
     raFunc.resize(3);
     decFunc.resize(3);
     getBestFitQuadratic(epoch, raFunc[0], raFunc[1], raFunc[2],
			 decFunc[0], decFunc[1], decFunc[2]);
     
}



void Track::calculateBestFitQuadratic(const std::vector<MopsDetection> &allDets, 
	       const bool useFullRaFit)
{
    int trackLen = componentDetectionIndices.size();

// A is a matrix with one row per MopsDetection, with the values of the fitting
// functions at the time of that detection.

    int raFuncLen;
    if (useFullRaFit && trackLen >= 7) {
	 raFuncLen = 5;
    } else {
	 raFuncLen = 3;
    }

    Eigen::MatrixXd raA(trackLen, raFuncLen);
    Eigen::VectorXd raCorr(trackLen);
    Eigen::MatrixXd decA(trackLen, 3);

// b is a vector with the measured values, either ra or dec, for each MopsDetection

    Eigen::VectorXd raB(trackLen);
    Eigen::VectorXd decB(trackLen);

// vectors of measurement errors

    Eigen::VectorXd raE(trackLen);
    Eigen::VectorXd decE(trackLen);


    int i = 0;
    for (std::set<unsigned int>::const_iterator detIndIt = componentDetectionIndices.begin();
         detIndIt != componentDetectionIndices.end(); detIndIt++, i++) {
        const MopsDetection* curDet = &allDets.at(*detIndIt);
	double t = curDet->getEpochMJD();
	double ra = curDet->getRA();
	double dec = curDet->getDec();
	double raTopoCorr = curDet->getRaTopoCorr();
	double raErr = curDet->getRaErr();
	double decErr = curDet->getDecErr();
	
	raA(i, 0) = 1.0;
	raA(i, 1) = t;
	if (raFuncLen==5) {
	     raCorr(i) = raTopoCorr;
	}
	raB(i) = ra;

	decA(i, 0) = 1.0;
	decB(i) = dec;

	raE(i) = raErr;
	decE(i) = decErr;

    }

// demean t, to reduce condition number.  Member 'epoch' is mean(t) - should rename

    epoch = raA.col(1).mean();
    raA.col(1).array() -= epoch;
    decA.col(1) = raA.col(1);

// calculate needed powers of t

    raA.col(2).array() = raA.col(1).array() * raA.col(1).array();
    decA.col(2) = raA.col(2);

// demean the raCorr column, if used, to avoid degeneracy with the constant term in the fit

    if (raFuncLen==5) {
	raA.col(3).array() = raA.col(2).array() * raA.col(1).array(); 
	raA.col(4).array() = raCorr.array() - raCorr.array().mean();
    }

// Solve in a least squares sense, raA * raFunc = raB, and similarly for dec
// raX and decX should be members of Track
// NOTE:  need to make this be weighted lsq - leave for the moment

    raFunc = raA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(raB);
    decFunc = decA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(decB);

// Calculate the residuals

    Eigen::VectorXd raResid = raB - raA * raFunc;
    Eigen::VectorXd raWtResid = raResid.array() / raE.array();
    chisqRa = raWtResid.dot(raWtResid);

    Eigen::VectorXd decResid = decB - decA * decFunc;
    Eigen::VectorXd decWtResid = decResid.array() / decE.array();
    chisqDec = decWtResid.dot(decWtResid);

// Calculate the prob(chisq), which will be the quality measure of the fit

    probChisqRa = gsl_cdf_chisq_Q(chisqRa, trackLen);
    probChisqDec = gsl_cdf_chisq_Q(chisqDec, trackLen);

#ifdef DEBUG
	      std::cout << "raB: \n" << raB << '\n';
	      std::cout << "raE: \n" << raE << '\n';
	      std::cout << "raA: \n" << raA << '\n';
	      std::cout << "raSVD: \n" << raA.jacobiSvd().singularValues() << '\n';
	      std::cout << "raFunc: \n" << raFunc << '\n';
	      std::cout << "raResid: \n" << raResid << '\n';
	      std::cout << "ra: chisq prob dof " << chisqRa << " " << probChisqRa << " " << trackLen << '\n';
	      std::cout << "decB: \n" << decB << '\n';
	      std::cout << "decE: \n" << decE << '\n';
	      std::cout << "decA: \n" << decA << '\n';
	      std::cout << "decSVD: \n" << decA.jacobiSvd().singularValues() << '\n';
	      std::cout << "decFunc: \n" << decFunc << '\n';
	      std::cout << "decResid: \n" << decResid << '\n';
	      std::cout << "dec: chisq prob dof " << chisqDec << " " << probChisqDec << " " << trackLen << '\n';
#endif
}



void Track::setNumUniqueNights(unsigned int n)
{
     numUniqueNights = n;
}
unsigned int Track::getNumUniqueNights() const
{
     return numUniqueNights;
}


void Track::predictLocationAtTime(const double mjd, double &ra, double &dec) const
{
/* 
   Note use of first guess at ra, used for calculating topocentric correction if needed.
   Then, final ra uses the whole formalism.
*/

    double t = mjd - epoch;

    Eigen::Vector3d tPowers(1.0, t, t*t);
    ra = raFunc.head(3).dot(tPowers);
    dec = decFunc.dot(tPowers);

    if (raFunc.size() == 5) {
	 MopsDetection tmpDet(0, mjd, ra, dec);
	 tmpDet.calculateTopoCorr();  
	 double raTopoCorr = tmpDet.getRaTopoCorr();
	 Eigen::Vector4d tPowersCubic(1.0, t, t*t, t*t*t);
	 ra = raFunc.head(4).dot(tPowersCubic) + raFunc(4)*raTopoCorr;
    }

#ifdef DEBUG
    std::cout << "tPowers: \n" << tPowers << '\n';
    std::cout << "raFunc: \n" << raFunc << '\n';
    std::cout << "ra: " << ra << '\n';
    std::cout << "decFunc: \n" << decFunc << '\n';
    std::cout << "dec: " << dec << '\n';
#endif
}



void Track::getBestFitQuadratic(double &epoch, double &ra0, double &raV, double &raAcc,
                                double &dec0, double &decV, double &decAcc) const
{

    epoch = this->epoch;
    ra0 = raFunc(0);
    raV = raFunc(1);
    raAcc = 2.0*raFunc(2);
    dec0 = decFunc(0);
    decV = decFunc(1);
    decAcc = 2.0*decFunc(2);

}

double Track::getFitRange() const {
     if (raFunc.size() == 5) {
	  return raFunc(4);
	  }
     else {
	  return 0;
     }
}


unsigned int Track::getId() const
{
     return myId;
}

void Track::setId(unsigned int newId)
{
     myId = newId;
}



MopsDetection Track::getLastDetection(
    const std::vector<MopsDetection> &dets)
    const
{
     if (componentDetectionIndices.size() < 1) {
	  LSST_EXCEPT(UninitializedException, 
   "Track: start time requested, but track contains no detections.");
     }
          
     bool haveLast = false;
     MopsDetection lastDet;
     double lastTime;

     std::set<unsigned int>::const_iterator indIt;
     for (indIt = componentDetectionIndices.begin(); 
	  indIt != componentDetectionIndices.end(); indIt++) {

	  if (*indIt > dets.size()) {
	       LSST_EXCEPT(BadIndexException, 
  "Track: the detection vector we were handed cannot possibly go with this track.");
	  }

	  double thisTime = dets.at(*indIt).getEpochMJD();
	  if ((!haveLast) || (thisTime > lastTime)) {
	       lastTime = thisTime;
	       haveLast = true;
               lastDet = dets.at(*indIt);
	  }
     }
     return lastDet;    
}




MopsDetection Track::getFirstDetection(
    const std::vector<MopsDetection> &dets) const
{
     if (componentDetectionIndices.size() < 1) {
	  LSST_EXCEPT(UninitializedException, 
   "Track: start time requested, but tracklet contains no detections.");
     }
          
     bool haveFirst = false;
     MopsDetection firstDet;
     double firstTime;

     std::set<unsigned int>::const_iterator indIt;
     for (indIt = componentDetectionIndices.begin(); 
	  indIt != componentDetectionIndices.end(); indIt++) {

	  if (*indIt > dets.size()) {
	       LSST_EXCEPT(BadIndexException, 
  "Track: the detection vector we were handed cannot possibly go with this track.");
	  }

	  double thisTime = dets.at(*indIt).getEpochMJD();
	  if ((!haveFirst) || (thisTime < firstTime)) {
	       firstTime = thisTime;
	       haveFirst = true;
               firstDet = dets.at(*indIt);
	  }
     }
     return firstDet;    
}



double Track::getStartTime(const std::vector<MopsDetection> &dets) const
{
    return getFirstDetection(dets).getEpochMJD();
}

double Track::getDeltaTime(const std::vector<MopsDetection> &dets) const
{
    return getLastDetection(dets).getEpochMJD() - getStartTime(dets);
}



}};

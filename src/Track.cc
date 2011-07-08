#include <algorithm>
#include <set>

#include "lsst/mops/Track.h"
#include "lsst/mops/Exceptions.h"

#undef DEBUG

namespace lsst { namespace mops {

Track::Track()
{
    epoch = 0;
}


void Track::addTracklet(unsigned int trackletIndex, 
			const Tracklet &t, 
			const std::vector<MopsDetection> & allDets)
{
     componentTrackletIndices.insert(trackletIndex);
     std::set<unsigned int>::const_iterator trackletDetIter;
     for (trackletDetIter = t.indices.begin(); 
	  trackletDetIter != t.indices.end(); 
	  trackletDetIter++) {
	  componentDetectionIndices.insert(*trackletDetIter);
	  componentDetectionDiaIds.insert(allDets.at(*trackletDetIter).getID());
     }

}
	  


void Track::addDetection(unsigned int detIndex, const std::vector<MopsDetection> & allDets)
{
    // TDB: might be nice to check that we haven't been given two dets from the same image time.
    // this will enforce some safety with quadratic fitting. not critical for now, 
    // since linkTracklets never does anything that dumb.
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
    epoch = other.epoch;
    raFunc = other.raFunc;
    decFunc = other.decFunc;
    return *this;
}




int Track::getObjectId(std::vector<MopsDetection> allDets) 
{
     if (componentDetectionIndices.size() == 0) {
	  throw LSST_EXCEPT(UninitializedException, 
			    "Cannot request object ID for track with no detections.");
     }
     // no object should have ssmId -1; it's -1 for noise >=0 for real ssmIds.
     int toRet = -2;
     std::set<unsigned int>::const_iterator detInd;
     for (detInd = componentDetectionIndices.begin();
	  detInd != componentDetectionIndices.end(); detInd++) {
	  int curId = allDets.at(*detInd).getSsmId();
	  if (toRet == -2) {
	       toRet = curId;
	  }
	  if ((curId == -1) || (toRet != curId))  {
	       return -1;
	  }	  
     }
     return toRet;
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


}};

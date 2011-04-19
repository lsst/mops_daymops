#include <algorithm>
#include <set>

#include "lsst/mops/Track.h"
#include "lsst/mops/Exceptions.h"

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








void Track::calculateBestFitQuadratic(const std::vector<MopsDetection> &allDets)
{
    int trackLen = componentDetectionIndices.size();

// A is a matrix with one row per MopsDetection, with the values of the fitting
// functions at the time of that detection.

    int raFuncLen;
    if (trackLen >= 5) {
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

    if (raFuncLen==5) {
	      std::cerr << "raB: \n" << raB << '\n';
	      std::cerr << "raE: \n" << raE << '\n';
	      std::cerr << "raA: \n" << raA << '\n';
	      std::cerr << "raSVD: \n" << raA.jacobiSvd().singularValues() << '\n';
	      std::cerr << "raFunc: \n" << raFunc << '\n';
	      std::cerr << "raResid: \n" << raResid << '\n';
	      std::cerr << "ra: chisq prob dof " << chisqRa << " " << probChisqRa << " " << trackLen << '\n';
	      std::cerr << "decB: \n" << decB << '\n';
	      std::cerr << "decE: \n" << decE << '\n';
	      std::cerr << "decA: \n" << decA << '\n';
	      std::cerr << "decSVD: \n" << decA.jacobiSvd().singularValues() << '\n';
	      std::cerr << "decFunc: \n" << decFunc << '\n';
	      std::cerr << "decResid: \n" << decResid << '\n';
	      std::cerr << "dec: chisq prob dof " << chisqDec << " " << probChisqDec << " " << trackLen << '\n';
	 }
}



void Track::predictLocationAtTime(const double mjd, double &ra, double &dec) const
{
/* 
   Want to avoid recalculating topocentric correction.  Can't just ignore it, unless 
   we're willing to run with much larger trackAdditionThreshold.  So, need to interpolate 
   it.
*/

// could make a MopsDetection with input epoch, calculate an approx ra/dec, call CalcTopoCorr()

// JUST QUADRATIC for the moment


    double t = mjd - epoch;

    Eigen::Vector3d tPowers(1.0, t, t*t);
    ra = raFunc.head(3).dot(tPowers);
    dec = decFunc.dot(tPowers);
#ifdef DEBUG
    std::cerr << "tPowers: \n" << tPowers << '\n';
    std::cerr << "raFunc: \n" << raFunc << '\n';
    std::cerr << "ra: " << ra << '\n';
    std::cerr << "decFunc: \n" << decFunc << '\n';
    std::cerr << "dec: " << dec << '\n';
#endif
}



void Track::getBestFitQuadratic(double &epoch, double &ra0, double &raV, double &raAcc,
                                double &dec0, double &decV, double &decAcc) const
{

    epoch = this->epoch;
    ra0 = raFunc(0);
    raV = raFunc(1);
    raAcc = raFunc(2);
    dec0 = decFunc(0);
    decV = decFunc(1);
    decAcc = decFunc(2);

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

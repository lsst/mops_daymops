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

    Eigen::MatrixXf raA(trackLen, 5);
    Eigen::MatrixXf decA(trackLen, 3);

// b is a vector with the measured values, either ra or dec, for each MopsDetection

    Eigen::VectorXf raB(trackLen);
    Eigen::VectorXf decB(trackLen);

// vectors of measurement errors

    Eigen::VectorXf raE(trackLen);
    Eigen::VectorXf decE(trackLen);


    int i = 0;
    for (std::set<unsigned int>::const_iterator detIndIt = componentDetectionIndices.begin();
         detIndIt != componentDetectionIndices.end(); detIndIt++, i++) {
        const MopsDetection* curDet = &allDets.at(*detIndIt);
	double t = curDet->getEpochMJD();
	double ra = curDet->getRA();
	double dec = curDet->getDec();
	double raCorr = curDet->getRaTopoCorr();
	double raErr = curDet->getRaErr();
	double decErr = curDet->getDecErr();
	
	double t2 = t*t;
	double t3 = t2*t;

	raA(i, 0) = 1.0;
	raA(i, 1) = t;
	raA(i, 2) = t2;
	raA(i, 3) = t3;
	raA(i, 4) = raCorr;
	raB(i) = ra;

	decA(i, 0) = 1.0;
	decA(i, 1) = t;
	decA(i, 2) = t2;
	decB(i) = dec;

	raE(i) = raErr;
	decE(i) = decErr;

    }

// demean the raCorr column, to avoid degeneracy with the constant term in the fit

    raA.col(4).array() -= raA.col(4).mean();

// Solve in a least squares sense, raA * raFunc = raB, and similarly for dec
// raX and decX should be members of Track
// NOTE:  need to make this be weighted lsq - leave for the moment

    raFunc = raA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(raB);
    decFunc = decA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(decB);

// Calculate the residuals

    Eigen::VectorXf raResid = raB - raA * raFunc;
    Eigen::VectorXf raWtResid = raResid.array() / raE.array();
    chisqRa = raWtResid.dot(raWtResid);

    Eigen::VectorXf decResid = decB - decA * decFunc;
    Eigen::VectorXf decWtResid = decResid.array() / decE.array();
    chisqDec = decWtResid.dot(decWtResid);

// Calculate the prob(chisq), which will be the quality measure of the fit
// NEED to be more careful with # of degrees of freedom

    probChisqRa = gsl_cdf_chisq_P(chisqRa, trackLen - 4);
    probChisqDec = gsl_cdf_chisq_P(chisqDec, trackLen - 2);

    std::cerr << "raB: \n" << raB << '\n';
    std::cerr << "raE: \n" << raE << '\n';
    std::cerr << "raA: \n" << raA << '\n';
    std::cerr << "raFunc: \n" << raFunc << '\n';
    std::cerr << "raResid: \n" << raResid << '\n';
    std::cerr << "ra: chisq prob dof " << chisqRa << " " << probChisqRa << " " << trackLen - 4 << '\n';

}



void Track::predictLocationAtTime(const double mjd, double &ra, double &dec) const
{
/*
    double t = mjd - epoch;
    ra = raFunc.at(0) + raFunc.at(1) * t + .5 * raFunc.at(2) * t * t;
    dec = decFunc.at(0) + decFunc.at(1) * t + .5 * decFunc.at(2) * t * t;
*/
}



void Track::getBestFitQuadratic(double &epoch, double &ra0, double &raV, double &raAcc,
                                double &dec0, double &decV, double &decAcc) const
{
/*
    epoch = this->epoch;
    ra0 = raFunc.at(0);
    raV = raFunc.at(1);
    raAcc = raFunc.at(2);
    dec0 = decFunc.at(0);
    decV = decFunc.at(1);
    decAcc = decFunc.at(2);
*/
}




}};

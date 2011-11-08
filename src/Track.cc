#include <algorithm>
#include <set>

#include "lsst/mops/Track.h"
#include "lsst/mops/Exceptions.h"

#define DEBUG

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
				      const bool useFullRaFit, std::ostream *outFile)
{
    int trackLen = componentDetectionIndices.size();

// A is a matrix with one row per MopsDetection, with the values of the fitting
// functions at the time of that detection.

    int raFuncLen, decFuncLen;
    if (useFullRaFit && trackLen >= 6) {
	 raFuncLen = 5;
	 decFuncLen = 4;
    } else {
	 raFuncLen = 3;
	 decFuncLen = 3;
    }

    Eigen::MatrixXd raA(trackLen, raFuncLen);
    Eigen::VectorXd raCorr(trackLen);
    Eigen::MatrixXd decA(trackLen, decFuncLen);

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
	decA.col(3) = raA.col(3);
    }

// Solve in a least squares sense, raA * raFunc = raB, and similarly for dec
// raX and decX should be members of Track
// NOTE:  need to make this be weighted lsq - leave for the moment

    Eigen::JacobiSVD<Eigen::MatrixXd> raSvd = raA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd raSingularValues = raSvd.singularValues();
    raFunc = raSvd.solve(raB);

    Eigen::JacobiSVD<Eigen::MatrixXd> decSvd = decA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd decSingularValues = decSvd.singularValues();
    decFunc = decSvd.solve(decB);

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

// Calculate the covariance matrices of the least squares solutions for ra and dec - based on Numerical
// Recipes eqn 15.4.20

    raCov.resize(raFuncLen, raFuncLen);
    Eigen::MatrixXd raV = raSvd.matrixV();
    for (int i=0; i<raFuncLen; i++) {
      for (int j=0; j<=i; j++) {
	double sum = 0;
	for (int k=0; k<raFuncLen; k++) {
	  sum += raV(i,k)*raV(j,k)/(raSingularValues(k)*raSingularValues(k));
	}
	raCov(i,j) = sum;
	if (j<i) raCov(j,i) = raCov(i,j);
      }
    }

    decCov.resize(decFuncLen, decFuncLen);
    Eigen::MatrixXd decV = decSvd.matrixV();
    for (int i=0; i<decFuncLen; i++) {
      for (int j=0; j<=i; j++) {
	double sum = 0;
	for (int k=0; k<decFuncLen; k++) {
	  sum += decV(i,k)*decV(j,k)/(decSingularValues(k)*decSingularValues(k));
	}
	decCov(i,j) = sum;
	if (j<i) decCov(j,i) = decCov(i,j);
      }
    }

#ifdef DEBUG

    if (outFile) {
	 double raCondNumber = raSingularValues(0)/raSingularValues(raSingularValues.rows()-1);
	 double decCondNumber = decSingularValues(0)/decSingularValues(decSingularValues.rows()-1);

	      *outFile << "raB: \n" << raB << '\n';
	      *outFile << "raE: \n" << raE << '\n';
	      *outFile << "raA: \n" << raA << '\n';
	      *outFile << "raSVD: \n" << raSingularValues << '\n';
	      *outFile << "raCov: \n" << raCov << '\n';
	      *outFile << "raFunc: \n" << raFunc << '\n';
	      *outFile << "raResid: \n" << raResid << '\n';
	      *outFile << "ra: chisq prob dof " << chisqRa << " " << probChisqRa << " " << trackLen << '\n';
	      *outFile << "ra npts, probChisq, condNum, : " << raB.rows() << " " << probChisqRa << " " << raCondNumber << '\n';
	      *outFile << "decB: \n" << decB << '\n';
	      *outFile << "decE: \n" << decE << '\n';
	      *outFile << "decA: \n" << decA << '\n';
	      *outFile << "decSVD: \n" << decSingularValues << '\n';
	      *outFile << "decCov: \n" << decCov << '\n';
	      *outFile << "decFunc: \n" << decFunc << '\n';
	      *outFile << "decResid: \n" << decResid << '\n';
	      *outFile << "dec: chisq prob dof " << chisqDec << " " << probChisqDec << " " << trackLen << '\n';
	      *outFile << "dec npts, probChisq, condNum, : " << decB.rows() << " " << probChisqDec << " " << decCondNumber << '\n';
    }
#endif
}



void Track::predictLocationAtTime(const double mjd, double &ra, double &dec) const
{
/* 
   Note use of first guess at ra, used for calculating topocentric correction if needed.
   Then, final ra uses the whole formalism.
*/

    double t = mjd - epoch;

    Eigen::Vector4d tPowersCubic(1.0, t, t*t, t*t*t);

    if (raFunc.size() == 5) {
         // initial guess at ra, dec to get topoCorr
         ra = raFunc.head(3).dot(tPowersCubic.head(3));
	 dec = decFunc.head(3).dot(tPowersCubic.head(3));
	 MopsDetection tmpDet(0, mjd, ra, dec);
	 tmpDet.calculateTopoCorr();  
	 double raTopoCorr = tmpDet.getRaTopoCorr();
	 ra = raFunc.head(4).dot(tPowersCubic) + raFunc(4)*raTopoCorr;
    } else {
	 ra = raFunc.head(3).dot(tPowersCubic.head(3));
    }

    if (decFunc.size() == 4) {
	 dec = decFunc.head(4).dot(tPowersCubic);
    } else {
	 dec = decFunc.head(3).dot(tPowersCubic.head(3));
    }

#ifdef DEBUG
    std::cout << "tPowers: \n" << tPowersCubic << '\n';
    std::cout << "raFunc: \n" << raFunc << '\n';
    std::cout << "ra: " << ra << '\n';
    std::cout << "decFunc: \n" << decFunc << '\n';
    std::cout << "dec: " << dec << '\n';
#endif
}

void Track::predictLocationUncertaintyAtTime(const double mjd, double &raUnc, double &decUnc) const
{

    double t = mjd - epoch;
    Eigen::VectorXd gVecRa;
    Eigen::VectorXd gVecDec;
    Eigen::MatrixXd tmp;

    if (raFunc.size() == 5) {
         double ra, dec;
	 predictLocationAtTime(mjd, ra, dec);
	 MopsDetection tmpDet(0, mjd, ra, dec);
	 tmpDet.calculateTopoCorr();  
	 double raTopoCorr = tmpDet.getRaTopoCorr();
	 gVecRa.resize(5);
	 gVecRa(0)=1.0;
	 gVecRa(1)=t;
	 gVecRa(2)=t*t;
	 gVecRa(3)=t*t*t;
	 gVecRa(4)=raTopoCorr;

    } else {
	 gVecRa.resize(3);
	 gVecRa(0)=1.0;
	 gVecRa(1)=t;
	 gVecRa(2)=t*t;
    }
    
    tmp = gVecRa * raCov * gVecRa.transpose();
    raUnc = tmp(0,0);

    
    if (decFunc.size() == 4) {
	 gVecDec.resize(4);
	 gVecDec(0)=1.0;
	 gVecDec(1)=t;
	 gVecDec(2)=t*t;
	 gVecDec(3)=t*t*t;

    } else {
	 gVecDec.resize(3);
	 gVecDec(0)=1.0;
	 gVecDec(1)=t;
	 gVecDec(2)=t*t;
    }
    tmp = gVecRa * raCov * gVecRa.transpose();
    decUnc = tmp(0,0);
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

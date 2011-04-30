// -*- LSST-C++ -*-


/* jmyers 2/11/11
 */

#include <algorithm>
#include <gsl/gsl_fit.h>

#include "lsst/mops/rmsLineFit.h"
#include "lsst/mops/Tracklet.h"
#include "lsst/mops/rmsLineFit.h"
#include "lsst/mops/common.h"

namespace lsst {
namespace mops {

Tracklet::Tracklet()
{
    myId = -1;
    isCollapsed = false; 
    RaBestFitFunc.clear(); 
    DecBestFitFunc.clear();
}



/* the final parameter is modified; it will hold Detections associated
   with the tracklet t. */

void Tracklet::getAllDetections(
    const std::vector<MopsDetection> & allDetections,
    std::vector<MopsDetection> &detectionsForTracklet) const
{
    detectionsForTracklet.clear();
    std::set<unsigned int>::const_iterator trackletIndexIter;

    for (trackletIndexIter = indices.begin();
         trackletIndexIter != indices.end();
         trackletIndexIter++) {
        detectionsForTracklet.push_back(allDetections.at(
                                            *trackletIndexIter));
    }
}







void Tracklet::calculateBestFitFunc(const std::vector<MopsDetection> &allDets) 
{
    /*
      jmyers: 

      TBD: this code is basically copy-pasted from rmsLineFit's
      leastSquaresSolveForRADecLinear, except that we calculate our
      own epoch and they use one specified by an outside user.  

      I did this initially because I was having SCons-related issues
      trying to get the darn thing to compile, ugh! But really, THIS
      is where this code belongs, so go back later and kill off the
      old function in rmsLineFit and replace all its calls with member
      calls to Tracklet.
     */

    std::vector<MopsDetection> trackletDets;
    getAllDetections(allDets, 
                     trackletDets);
    
    unsigned int numDets = trackletDets.size();
    RaBestFitFunc.clear();
    DecBestFitFunc.clear();

    fitFuncEpoch = getStartTime(allDets);
    
    std::vector<double> MJDs(numDets);
    std::vector<double> RAs(numDets);
    std::vector<double> Decs(numDets);
    
    
    for (unsigned int i = 0; i < numDets; i++) {
        RAs[i] = trackletDets[i].getRA();
        Decs[i] = trackletDets[i].getDec();
        MJDs[i] = trackletDets[i].getEpochMJD() - fitFuncEpoch;
    }
    
    if (*std::max_element(RAs.begin(), RAs.end()) - 
        *std::min_element(RAs.begin(), RAs.end()) > 180.) {
        make180To360Negative(RAs);
    }
    if (*std::max_element(Decs.begin(), Decs.end()) - 
        *std::min_element(Decs.begin(), Decs.end()) > 180.) {
        make180To360Negative(Decs);
    }
    if ((*std::max_element(RAs.begin(), RAs.end()) - 
         *std::min_element(RAs.begin(), RAs.end()) > 180.) ||
        (*std::max_element(Decs.begin(), Decs.end()) - 
         *std::min_element(Decs.begin(), Decs.end()) > 180.)) {
        throw LSST_EXCEPT(ProgrammerErrorException, 
"EE: Unexpected coding error: could not move data into a contiguous < 180 degree range.\n");
    }
    
    /* use linear least squares to get MJD->RA, MJD->Dec funs */
    /* need to get C-style arrays for gsl. Do this The Right Way, rather
     * than with the famous hack: &(RA[0]) */
    double *arrayRAs = (double*)malloc(sizeof(double) * numDets);
    double *arrayDecs = (double*)malloc(sizeof(double) * numDets);
    double *arrayMJDs = (double*)malloc(sizeof(double) * numDets);
    
    if ((arrayRAs == NULL) || (arrayDecs == NULL) || (arrayMJDs == NULL)) {
        throw LSST_EXCEPT(pexExcept::MemoryException, 
"Malloc returned NULL on a very small malloc. System out of memory or something very odd.\n");
    }
    
    for (unsigned int i = 0; i < numDets; i++) {
        arrayRAs[i] = RAs[i];
        arrayDecs[i] = Decs[i];
        arrayMJDs[i] = MJDs[i];
    }
    
    int gslRV = 0;
    double offset, slope;
    double cov00, cov01, cov11, sumsq;
    
    // arrayMJDs is independent, arrayRAs is dependent.
    gslRV += gsl_fit_linear(arrayMJDs, 1, arrayRAs, 1, numDets, &offset, &slope, 
                            &cov00,&cov01,&cov11,&sumsq);
    RaBestFitFunc.push_back(offset);
    RaBestFitFunc.push_back(slope);
    
    gslRV += gsl_fit_linear(arrayMJDs, 1, arrayDecs, 1, numDets, &offset, &slope, 
                            &cov00,&cov01,&cov11,&sumsq);
    DecBestFitFunc.push_back(offset);
    DecBestFitFunc.push_back(slope);
    
    if (gslRV != 0) {
        throw LSST_EXCEPT(GSLException, 
                          "EE: gsl_fit_linear unexpectedly returned error.\n");
    }
    
    free(arrayRAs);
    free(arrayDecs);
    free(arrayMJDs);
}




Tracklet::Tracklet(std::set <unsigned int> startIndices) 
{ 
    myId = -1;
    RaBestFitFunc.clear();
    DecBestFitFunc.clear();
    isCollapsed = false; 
    indices = startIndices;
}
     







void Tracklet::getFitFunction(double& epoch, std::vector<double> &raFunc,
                              std::vector<double> &decFunc) 
{
    if ((RaBestFitFunc.size() == 0) || (DecBestFitFunc.size() == 0)) {
        throw LSST_EXCEPT(UninitializedException,
" Tracklet got request for best-fit function, but never got request to calculate one!");
    }
    
    epoch = fitFuncEpoch;
    raFunc = RaBestFitFunc;
    decFunc = DecBestFitFunc;

}






MopsDetection Tracklet::getFirstDetection(
    const std::vector<MopsDetection> &dets) const
{
     if (indices.size() < 1) {
	  LSST_EXCEPT(UninitializedException, 
   "Tracklet: start time requested, but tracklet contains no detections.");
     }
          
     bool haveFirst = false;
     MopsDetection firstDet;
     double firstTime;

     std::set<unsigned int>::const_iterator indIt;
     for (indIt = indices.begin(); indIt != indices.end(); indIt++) {

	  if (*indIt > dets.size()) {
	       LSST_EXCEPT(BadIndexException, 
  "Tracklet: the detection vector we were handed cannot possibly go with this tracklet.");
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



double Tracklet::getStartTime(const std::vector<MopsDetection> &dets) const
{
    return getFirstDetection(dets).getEpochMJD();
}



double Tracklet::getDeltaTime(const std::vector<MopsDetection> &dets) const
{
    return getLastDetection(dets).getEpochMJD() - getStartTime(dets);
}


MopsDetection Tracklet::getLastDetection(
    const std::vector<MopsDetection> &dets)
    const
{
     if (indices.size() < 1) {
	  LSST_EXCEPT(UninitializedException, 
   "Tracklet: start time requested, but tracklet contains no detections.");
     }
          
     bool haveLast = false;
     MopsDetection lastDet;
     double lastTime;

     std::set<unsigned int>::const_iterator indIt;
     for (indIt = indices.begin(); indIt != indices.end(); indIt++) {

	  if (*indIt > dets.size()) {
	       LSST_EXCEPT(BadIndexException, 
  "Tracklet: the detection vector we were handed cannot possibly go with this tracklet.");
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


Tracklet & Tracklet::operator= (const Tracklet &other)
 {
    copy(other);
    return *this;
 } 

bool Tracklet::operator==(const Tracklet &other) const
{
    if (indices == other.indices){
        return true;
    }
    else {
        return false;
    }
}

bool Tracklet::operator!=(const Tracklet &other) const
{
    if (indices != other.indices){
        return true;
    }
    else {
        return false;
    }
    
}

bool Tracklet::operator<(const Tracklet &other) const
{
    return (indices < other.indices);
}

void Tracklet::copy(const Tracklet &other)
{
    indices = other.indices;
    RaBestFitFunc = other.RaBestFitFunc;
    DecBestFitFunc = other.DecBestFitFunc;
    isCollapsed = other.isCollapsed;
    myId = other.myId;
}

}}

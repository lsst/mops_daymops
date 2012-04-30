// -*- LSST-C++ -*-


/* jmyers 8/18/08
 *
 *This code is for getting the best-fit line for a series of observations
 * (i.e. a tracklet) and determining the RMS of the error for the detections.
 *
 * This also includes code for filtering out groups of tracklets which do not
 * fit the criteria.
 */


#ifndef LSST_RMSLINEFIT_H
#define LSST_RMSLINEFIT_H

#include <vector>
#include <map>

#include "lsst/mops/Tracklet.h"
#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/Exceptions.h"

namespace lsst {
namespace mops {


    /*
     * return the root mean squared distance of the tracklet from the line.
     * (i.e. the sqrt of the sum of the squared distances of each detection from
     * the best-fit line, all divided by the size of the tracklet).
     *
     * if perDetSqDist != NULL (the default value), then we ASSUME that
     * perDetSqDist is a pointer to an allocated vector, and we will populate it
     * s.t. mean Sq. distance (trackletDets[i],line) = perDetSqDist[i] 
     */
    double rmsForTracklet(Tracklet t, const std::vector<MopsDetection> *detections, 
                          std::vector<double>*perDetSqDist=NULL);
    


    
    /*
      given an array of detections, solve for RA and Dec in terms of MJD on the
      assumption that both are linear.  Output vectors will look like this when done:
      
      if RA = MJD * m + b
      RASlopeAndOffsetOut = [m, b] 
      
      (and similarly for Dec.)
      
      if timeOffset is specified as non-zero, this value will be substracted from the 
      MJDs of each detection.
    */
    void leastSquaresSolveForRADecLinear(const std::vector <MopsDetection> *trackletDets,
                                         std::vector<double> &RASlopeAndOffsetOut,
                                         std::vector<double> &DecSlopeAndOffsetOut, 
                                         double timeOffset=0.0);
    
/*
 * given a vector of Tracklets and the corresponding vector of MopsDetections,
 * add to output only those tracklets for which rms < maxRMSm * av. magnitude + maxRMSm
 */
    
    void filterByLineFitAddToOutputVector(const std::vector<Tracklet> *tracklets, 
                                          const std::vector<MopsDetection> * allDets,
                                          double maxRMS,
                                          std::vector<Tracklet> &output);
    
    // fetch all the detections and put them in a fresh vector.
    std::vector<MopsDetection> getTrackletDets(const Tracklet *t, 
                                               const std::vector<MopsDetection>* allDets);

    /*
     * given a pair of functions relating MJD to RA, Dec s.t. RA =
     * RASlope*(MJD - tOffset) + RAintercept and the similarly for
     * declination, return a map which relates indices of *t to
     * squared distance to the projected point.
     */
    std::map <unsigned int, double> 
    getPerDetSqDistanceToLine(const Tracklet *t, const std::vector<MopsDetection>* allDets, 
                              double RASlope, double RAIntercept, 
                              double DecSlope, double DecIntercept, 
                              double tOffset);




}} // close lsst::mops

#endif

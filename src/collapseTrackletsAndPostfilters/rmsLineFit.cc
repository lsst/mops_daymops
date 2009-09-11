// -*- LSST-C++ -*-

#include <cmath>
#include <fstream>
#include <iostream>

#include "TrackletCollapser.h"
#include "../Exceptions.h"
#include "rmsLineFit.h"


namespace ctExcept = collapseTracklets::exceptions;

namespace rmsLineFit {





    double getAverageMagnitude(const Tracklet t, const::std::vector<Detection>* detections) {
        double sum = 0.;
        unsigned int count = 0;
        std::set<unsigned int>::const_iterator detIndexIter;
        for (detIndexIter = t.indices.begin(); detIndexIter != t.indices.end(); detIndexIter++) {
            sum += (*detections)[*detIndexIter].getMag();
            count++;
        }
        if (count < 1) {
            throw LSST_EXCEPT(ctExcept::BadParameterException, "EE: getAverageMagnitude: highly unexpected error - tracklet has no detections?\n");
        }
        return sum/count;
    }




    double rmsForTracklet(const Tracklet t, const std::vector<Detection> *detections, 
                          std::vector<double> *perDetSqDist) {
        std::vector<double> RASlopeAndOffset, DecSlopeAndOffset;
        std::set<unsigned int>::iterator indicesIter;
        std::vector<Detection> trackletDets;
        std::vector<Detection>::iterator detIter;
        for (indicesIter = t.indices.begin(); 
             indicesIter != t.indices.end();
             indicesIter++) {
            trackletDets.push_back((*detections)[*indicesIter]);
        }
        
        if ((perDetSqDist != NULL) && (perDetSqDist->size() != 0)) {
            throw LSST_EXCEPT(ctExcept::ProgrammerErrorException, 
                              "EE: PROGRAMMING ERROR: perDetSqDist argument to rmsForTracklet must be either NULL or a pointer to an allocated, empty vector!\n");
        }
        

        if (t.indices.size() == 0) {
            throw LSST_EXCEPT(ctExcept::ProgrammerErrorException, 
                              "EE: PROGRAMMING ERROR: rmsForTracklet received a tracklet associated with 0 detection indices.\n");
        }
        
        /* set the first detection time as t = 0 (timeOffset) */
        double t0 = trackletDets[0].getEpochMJD();
        collapseTracklets::TrackletCollapser myTC;
        myTC.leastSquaresSolveForRADecLinear(&trackletDets, RASlopeAndOffset, 
                                                           DecSlopeAndOffset, t0);
        /* since our current least squares library doesn't give it to us, we have to re-calculate
         * what that least square was...*/
        double squaresSum = 0.0;
        for (detIter = trackletDets.begin(); detIter != trackletDets.end(); detIter++) {
            double tOffset = detIter->getEpochMJD() - t0;

            double projectedRA = RASlopeAndOffset[0] * tOffset + RASlopeAndOffset[1];
	    projectedRA = KDTree::Common::convertToStandardDegrees(projectedRA);

            double projectedDec = DecSlopeAndOffset[0] * tOffset + DecSlopeAndOffset[1];
	    projectedDec = KDTree::Common::convertToStandardDegrees(projectedDec);

	    double actualRA = KDTree::Common::convertToStandardDegrees(detIter->getRA());
	    double actualDec = KDTree::Common::convertToStandardDegrees(detIter->getDec());

	    // 	    std::cout << " projected position: (" << projectedRA << ", " << projectedDec << ")" << std::endl;
	    // 	    std::cout << " Actual position: (" << actualRA << ", " << actualDec <<  ")" << std::endl;

            double RADist = fabs(projectedRA -  actualRA);
            double DecDist = fabs(projectedDec -  actualDec);
	    if (RADist > 180.) {
	      RADist = fabs(RADist - 360);
	    }
	    if (DecDist > 180.) {
	      DecDist = fabs(DecDist - 360);
	    }
	    if ((RADist > 180.)  || (DecDist > 180.)) {
                throw LSST_EXCEPT(ctExcept::ProgrammerErrorException,
                                  "EE: rmsLineFit: Unexpected programming error\n");
            }
	    double localDistSquared = RADist*RADist + DecDist*DecDist;
            /* get the sum of all the squares of distances */
            if (perDetSqDist != NULL) {
                perDetSqDist->push_back(localDistSquared);
            }
            squaresSum += localDistSquared;
        }
	double result = sqrt(squaresSum) / trackletDets.size(); 
	return  result;
    }



    


    
    void filterByLineFitAddToOutputVector(const std::vector<Tracklet> *tracklets, 
                                          const std::vector<Detection> *allDets,
                                          double maxRMSm, double maxRMSb,
                                          std::vector<Tracklet> &output) {
        std::vector<Tracklet>::const_iterator tIter;
        for (tIter = tracklets->begin(); tIter != tracklets->end(); tIter++) {
	  double rms = rmsForTracklet(*tIter, allDets);
	  // std::cout << "got RMS = " << rms << ", tracklet length = " <<  tIter->indices.size() << std::endl;
	  double maxRMS =  maxRMSm * getAverageMagnitude(*tIter, allDets) + maxRMSb;
	  if ( maxRMS >= rms) {
                output.push_back(*tIter);
            }
        }        
    }




    
    /*
     * given a pair of functions relating MJD to RA, Dec s.t. RA = RASlope*(MJD - tOffset) + RAintercept and
     * the similarly for declination, return a map which relates indices of *t to squared distance to the projected point.
     */
    std::map <unsigned int, double> getPerDetSqDistanceToLine(const Tracklet *t, const std::vector<Detection>* allDets, 
                                                              double RASlope, double RAIntercept, double DecSlope, double DecIntercept, 
                                                              double tOffset) {
        std::map <unsigned int, double> results;
        for (std::set<unsigned int>::const_iterator indicesIter = t->indices.begin();
             indicesIter != t->indices.end();
             indicesIter++) {
            double t = (*allDets)[*indicesIter].getEpochMJD() - tOffset;
            double projectedRA = KDTree::Common::convertToStandardDegrees(RASlope * t + RAIntercept);
            double projectedDec = KDTree::Common::convertToStandardDegrees(DecSlope *t + DecIntercept);
            double obsRA = KDTree::Common::convertToStandardDegrees((*allDets)[*indicesIter].getRA());
            double obsDec = KDTree::Common::convertToStandardDegrees((*allDets)[*indicesIter].getDec());
            double RADist = fabs(projectedRA - obsRA);
            if (RADist > 180.) {
                RADist = fabs(RADist - 360);
                if (RADist > 180.) {
                    throw LSST_EXCEPT(ctExcept::ProgrammerErrorException,
                                      "EE: getPerDetSqDistanceToLine: Unexpected programming error: failed to get distance < 360 between two degree points\n");
                }
            }
            double DecDist = fabs(projectedDec - obsDec);
            if (DecDist > 180.) {
                DecDist = fabs(DecDist - 360);
                if (DecDist > 180.) {
                    throw LSST_EXCEPT(ctExcept::ProgrammerErrorException,
                                      "EE: Unexpected programming error: failed to get distance < 360 between two degree points\n");
                }
            }
            double sqDist = RADist*RADist + DecDist*DecDist;
            results.insert(std::make_pair(*indicesIter, sqDist));            
        }
        return results;
    }



    std::vector<Detection> getTrackletDets(const Tracklet *t, const std::vector<Detection>* allDets) {
        std::vector<Detection> results;
        for (std::set<unsigned int>::const_iterator iIter = t->indices.begin();
             iIter != t->indices.end();
             iIter++) {
            results.push_back((*allDets)[*iIter]);
        }
        return results;
    }



    
    Tracklet TrackletPurifier::purifyTracklet(const Tracklet *t, const std::vector<Detection>* allDets, 
                                              double maxRMSm, double maxRMSb) {
        Tracklet curTracklet = *t;
        bool isClean = false;

        while (isClean == false) {
            double t0 = (*allDets)[*(curTracklet.indices.begin())].getEpochMJD();
            std::vector<Detection> curTrackletDets = getTrackletDets(&curTracklet, allDets);
            std::vector<double> RASlopeAndOffset, DecSlopeAndOffset;
            collapseTracklets::TrackletCollapser myTC;
            myTC.leastSquaresSolveForRADecLinear(&curTrackletDets, RASlopeAndOffset, 
                                                 DecSlopeAndOffset, t0);
            std::map<unsigned int, double> indexToSqDist = 
                getPerDetSqDistanceToLine(&curTracklet, allDets, RASlopeAndOffset[0], RASlopeAndOffset[1],
                                          DecSlopeAndOffset[0], DecSlopeAndOffset[1], t0);

            isClean = true;
            double worstDetVal = 0.0;
            unsigned int worstDetIndex = 0;
            for (std::map<unsigned int, double>::iterator distIter = indexToSqDist.begin();
                 distIter != indexToSqDist.end(); distIter++) {
                double distMax = maxRMSm * (*allDets)[distIter->first].getMag() + maxRMSb;
                if ((distIter->second > distMax*distMax) && (distIter->second > worstDetVal)) {
                    worstDetVal = distIter->second;
		    if (worstDetVal > 1) {
                        std::cerr << "Warning: detection point to projected point is improbably large distance: " << worstDetVal << std::endl;
		    }
                    worstDetIndex = distIter->first;
                    isClean = false;
                }
            }
            if (isClean == false) {
                curTracklet.indices.erase(worstDetIndex);
            }
        }
        return curTracklet;
    }


    void TrackletPurifier::purifyTracklets(const std::vector<Tracklet> *trackletsVector,
                                           const std::vector<Detection> *detsVector,
                                           double maxRMSm, double maxRMSb, unsigned int minObs,
                                           std::vector<Tracklet> &output)
    {
        if (output.size() != 0) {
            throw LSST_EXCEPT(ctExcept::BadParameterException, 
                              "purifyTracklets: output vector not empty\n");
        }
        
        std::vector<Tracklet>::const_iterator tIter;
        for (tIter = trackletsVector->begin(); tIter != trackletsVector->end(); tIter++) {
            Tracklet tmp = purifyTracklet(&(*tIter), detsVector, maxRMSm, maxRMSb);
            if (tmp.indices.size() >= minObs) {
                output.push_back(tmp);
            }
        }        
    
        
    }


}

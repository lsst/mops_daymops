// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>

#include <gsl/gsl_fit.h>

#include "lsst/mops/rmsLineFit.h"
#include "lsst/mops/KDTree.h"




namespace lsst {
    namespace mops {



    void make180To360Negative(std::vector<double> &v) {
        std::vector<double>::iterator iter;
        for (iter = v.begin(); iter != v.end(); iter++){
            if (*iter > 180) {
                *iter = *iter - 360;
            }
        }
    }




    void leastSquaresSolveForRADecLinear(const std::vector <MopsDetection> *trackletDets,
                                         std::vector<double> &RASlopeAndOffsetOut,
                                         std::vector<double> &DecSlopeAndOffsetOut,
                                         double timeOffset) {
        unsigned int numDets = (*trackletDets).size();
        
        std::vector<double> MJDs(numDets);
        std::vector<double> RAs(numDets);
        std::vector<double> Decs(numDets);

        if ((RASlopeAndOffsetOut.size() != 0) || (DecSlopeAndOffsetOut.size() != 0)) {
            throw LSST_EXCEPT(BadParameterException, 
                              "EE:  leastSquaresSolveForRADecLinear: output vectors were not empty.\n");
        }

        for (unsigned int i = 0; i < numDets; i++) {
            RAs[i] = (*trackletDets)[i].getRA();
            Decs[i] = (*trackletDets)[i].getDec();
            MJDs[i] = (*trackletDets)[i].getEpochMJD() - timeOffset;
        }

        /* here we make the assumption that the data is not very far apart in RA
         * or Dec degrees.  Since we are currently only concerned with a single
         * night of observation, this should be very safe.  Note that under this
         * assumption, our data could only span > 180 deg if the object crosses
         * the 0/360 line.*/

        if (*std::max_element(RAs.begin(), RAs.end()) - *std::min_element(RAs.begin(), RAs.end()) > 180.) {
            make180To360Negative(RAs);
        }
        if (*std::max_element(Decs.begin(), Decs.end()) - *std::min_element(Decs.begin(), Decs.end()) > 180.) {
            make180To360Negative(Decs);            
        }
        if ((*std::max_element(RAs.begin(), RAs.end()) - *std::min_element(RAs.begin(), RAs.end()) > 180.) ||
            (*std::max_element(Decs.begin(), Decs.end()) - *std::min_element(Decs.begin(), Decs.end()) > 180.)) {
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
            throw LSST_EXCEPT(pexExcept::MemoryException, "Malloc returned NULL on a very small malloc. System out of memory or something very odd.\n");
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
        RASlopeAndOffsetOut.push_back(slope);
        RASlopeAndOffsetOut.push_back(offset);
        
        gslRV += gsl_fit_linear(arrayMJDs, 1, arrayDecs, 1, numDets, &offset, &slope, 
                                &cov00,&cov01,&cov11,&sumsq);
        DecSlopeAndOffsetOut.push_back(slope);
        DecSlopeAndOffsetOut.push_back(offset);

        if (gslRV != 0) {
            throw LSST_EXCEPT(GSLException, "EE: gsl_fit_linear unexpectedly returned error.\n");
        }

        free(arrayRAs);
        free(arrayDecs);
        free(arrayMJDs);

    }







    double rmsForTracklet(const Tracklet t, const std::vector<MopsDetection> *detections, 
                          std::vector<double> *perDetSqDist) {
        std::vector<double> RASlopeAndOffset, DecSlopeAndOffset;
        std::set<unsigned int>::iterator indicesIter;
        std::vector<MopsDetection> trackletDets;
        std::vector<MopsDetection>::iterator detIter;
        for (indicesIter = t.indices.begin(); 
             indicesIter != t.indices.end();
             indicesIter++) {
            trackletDets.push_back((*detections)[*indicesIter]);
        }
        
        if ((perDetSqDist != NULL) && (perDetSqDist->size() != 0)) {
            throw LSST_EXCEPT(ProgrammerErrorException, 
                              "EE: PROGRAMMING ERROR: perDetSqDist argument to rmsForTracklet must be either NULL or a pointer to an allocated, empty vector!\n");
        }
        

        if (t.indices.size() == 0) {
            throw LSST_EXCEPT(ProgrammerErrorException, 
                              "EE: PROGRAMMING ERROR: rmsForTracklet received a tracklet associated with 0 detection indices.\n");
        }
        
        /* set the first detection time as t = 0 (timeOffset) */
        double t0 = trackletDets[0].getEpochMJD();
        leastSquaresSolveForRADecLinear(&trackletDets, RASlopeAndOffset, 
                                        DecSlopeAndOffset, t0);
        /* since our current least squares library doesn't give it to us, we have to re-calculate
         * what that least square was...*/
        double squaresSum = 0.0;
        for (detIter = trackletDets.begin(); detIter != trackletDets.end(); detIter++) {
            double tOffset = detIter->getEpochMJD() - t0;
            
            double projectedRA = RASlopeAndOffset[0] * tOffset + RASlopeAndOffset[1];
	    projectedRA = convertToStandardDegrees(projectedRA);
            
            double projectedDec = DecSlopeAndOffset[0] * tOffset + DecSlopeAndOffset[1];
	    projectedDec = convertToStandardDegrees(projectedDec);

	    double actualRA = convertToStandardDegrees(detIter->getRA());
	    double actualDec = convertToStandardDegrees(detIter->getDec());

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
                throw LSST_EXCEPT(ProgrammerErrorException,
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
                                          const std::vector<MopsDetection> *allDets,
                                          double maxRMS,
                                          std::vector<Tracklet> &output) {
        std::vector<Tracklet>::const_iterator tIter;
        for (tIter = tracklets->begin(); tIter != tracklets->end(); tIter++) {
	  double rms = rmsForTracklet(*tIter, allDets);
	  // std::cout << "got RMS = " << rms << ", tracklet length = " <<  tIter->indices.size() << std::endl;
	  if ( maxRMS >= rms) {
                output.push_back(*tIter);
            }
        }        
    }




    
    /*
     * given a pair of functions relating MJD to RA, Dec s.t. RA = RASlope*(MJD - tOffset) + RAintercept and
     * the similarly for declination, return a map which relates indices of *t to squared distance to the projected point.
     */
    std::map <unsigned int, double> getPerDetSqDistanceToLine(const Tracklet *t, const std::vector<MopsDetection>* allDets, 
                                                              double RASlope, double RAIntercept, double DecSlope, double DecIntercept, 
                                                              double tOffset) {
        std::map <unsigned int, double> results;
        for (std::set<unsigned int>::const_iterator indicesIter = t->indices.begin();
             indicesIter != t->indices.end();
             indicesIter++) {
            double t = (*allDets)[*indicesIter].getEpochMJD() - tOffset;
            double projectedRA = convertToStandardDegrees(RASlope * t + RAIntercept);
            double projectedDec = convertToStandardDegrees(DecSlope *t + DecIntercept);
            double obsRA = convertToStandardDegrees((*allDets)[*indicesIter].getRA());
            double obsDec = convertToStandardDegrees((*allDets)[*indicesIter].getDec());
            double RADist = fabs(projectedRA - obsRA);
            if (RADist > 180.) {
                RADist = fabs(RADist - 360);
                if (RADist > 180.) {
                    throw LSST_EXCEPT(ProgrammerErrorException,
                                      "EE: getPerDetSqDistanceToLine: Unexpected programming error: failed to get distance < 360 between two degree points\n");
                }
            }
            double DecDist = fabs(projectedDec - obsDec);
            if (DecDist > 180.) {
                DecDist = fabs(DecDist - 360);
                if (DecDist > 180.) {
                    throw LSST_EXCEPT(ProgrammerErrorException,
                                      "EE: Unexpected programming error: failed to get distance < 360 between two degree points\n");
                }
            }
            double sqDist = RADist*RADist + DecDist*DecDist;
            results.insert(std::make_pair(*indicesIter, sqDist));            
        }
        return results;
    }



    std::vector<MopsDetection> getTrackletDets(const Tracklet *t, const std::vector<MopsDetection>* allDets) {
        std::vector<MopsDetection> results;
        for (std::set<unsigned int>::const_iterator iIter = t->indices.begin();
             iIter != t->indices.end();
             iIter++) {
            results.push_back((*allDets)[*iIter]);
        }
        return results;
    }



    
    Tracklet TrackletPurifier::purifyTracklet(const Tracklet *t, const std::vector<MopsDetection>* allDets, 
                                              double maxRMS) {
        Tracklet curTracklet = *t;
        bool isClean = false;

        while (isClean == false) {
            double t0 = (*allDets)[*(curTracklet.indices.begin())].getEpochMJD();
            std::vector<MopsDetection> curTrackletDets = getTrackletDets(&curTracklet, allDets);
            std::vector<double> RASlopeAndOffset, DecSlopeAndOffset;
            leastSquaresSolveForRADecLinear(&curTrackletDets, RASlopeAndOffset, 
                                                 DecSlopeAndOffset, t0);
            std::map<unsigned int, double> indexToSqDist = 
                getPerDetSqDistanceToLine(&curTracklet, allDets, RASlopeAndOffset[0], RASlopeAndOffset[1],
                                          DecSlopeAndOffset[0], DecSlopeAndOffset[1], t0);

            isClean = true;
            double worstDetVal = 0.0;
            unsigned int worstDetIndex = 0;
            for (std::map<unsigned int, double>::iterator distIter = indexToSqDist.begin();
                 distIter != indexToSqDist.end(); distIter++) {
                double distMax = maxRMS;
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
                                           const std::vector<MopsDetection> *detsVector,
                                           double maxRMS, unsigned int minObs,
                                           std::vector<Tracklet> &output)
    {
        if (output.size() != 0) {
            throw LSST_EXCEPT(BadParameterException, 
                              "purifyTracklets: output vector not empty\n");
        }
        
        std::vector<Tracklet>::const_iterator tIter;
        for (tIter = trackletsVector->begin(); tIter != trackletsVector->end(); tIter++) {
            Tracklet tmp = purifyTracklet(&(*tIter), detsVector, maxRMS);
            if (tmp.indices.size() >= minObs) {
                output.push_back(tmp);
            }
        }        
    
        
    }


    }} // close lsst::mops namespace

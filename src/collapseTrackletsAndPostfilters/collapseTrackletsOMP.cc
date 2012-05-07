// -*- LSST-C++ -*-
/* jonathan myers */


#include <cmath>
#include <istream>
#include <sstream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "lsst/mops/daymops/collapseTrackletsAndPostfilters/collapseTracklets.h"
#include "lsst/mops/rmsLineFit.h"

// these just need to be greater or less than any date we'll ever see...
#define IMPOSSIBLY_EARLY_MJD -1.0E16
#define IMPOSSIBLY_LATE_MJD 1.0E16

// hard to guess an ideal number for this - could use some benchmarks.
#define MAX_LEAF_SIZE 50 


namespace lsst {
    namespace mops {

    // declare some helper functions - these don't need to be exposed.

    void populateTrackletsForTreeVector(const std::vector<MopsDetection> *detections,
                                        const std::vector<Tracklet> *tracklets,
                                        std::vector<PointAndValue <unsigned int> >
                                        &trackletsForTree);




    /* returns true iff two tracklets are 'compatible', i.e. the set of their
       combined unique detections contains no detections with the same
       observation time. NB: if t1 is detections AB, and t2 is detections AC,
       then they ARE compatible if B, C have different times (because their
       combined detections are ABC). */

    bool trackletsAreCompatible(const std::vector<MopsDetection> * detections,
                                const Tracklet &t1, const Tracklet &t2);



    /* returns the average of the squared distances of each detection in queryTracklet from
     * the line described by slopeAndOffsetRA, slopeAndOffsetDec. */ 
    double getAverageSqDist(std::vector<double> slopeAndOffsetRA, std::vector<double> slopeAndOffsetDec,
                            const std::vector<MopsDetection>* allDets, const Tracklet* queryTracklet, 
                            double timeOffset=0.0);



    void printTracklet(Tracklet t) {
        std::cout << "collapsed: " << t.isCollapsed << std::endl;
        for (std::set<unsigned int>::iterator iter = t.indices.begin();
             iter != t.indices.end(); iter++) {
            std::cout << *iter << " ";
        }
        std::cout << std::endl;
    }












    // definitions of our functions






    
    void printTrackletVector(std::vector<Tracklet> vt) {
        for (std::vector<Tracklet>::iterator iter = vt.begin();
             iter != vt.end(); iter++) {
            printTracklet(*iter);
        }
    }





    /* returns true iff smaller is a subset of larger*/
    bool isASubset(const std::set<unsigned int>* smaller, 
                    const std::set<unsigned int>* larger) {
        if (smaller->size() > larger->size()) {
            return false;
        }

        std::set<unsigned int> uniqueIndices;

        std::set_union(smaller->begin(), smaller->end(), 
                       larger->begin(), larger->end(),
                       std::inserter(uniqueIndices, uniqueIndices.begin()));
        if (uniqueIndices.size() == larger->size()) {
            return true;
        }
        return false;
        
    }



    double getSqDist(std::vector<double> slopeAndOffsetRA, std::vector<double> slopeAndOffsetDec,
                     const MopsDetection* det, double timeOffset) {
        double mjd = det->getEpochMJD();
        double projectedRA = slopeAndOffsetRA[0]*(mjd - timeOffset) + slopeAndOffsetRA[1];
        double projectedDec = slopeAndOffsetDec[0]*(mjd - timeOffset) + slopeAndOffsetDec[1];
        double obsRA = det->getRA();
        double obsDec = det->getDec();
        double distRA = obsRA - projectedRA;
        double distDec = obsDec - projectedDec;
        return (distRA)*(distRA) + (distDec)*(distDec);
        
    }




    double getAverageSqDist(std::vector<double> slopeAndOffsetRA, std::vector<double> slopeAndOffsetDec,
                            const std::vector<MopsDetection>* allDets, const Tracklet* queryTracklet, 
                            double timeOffset) {
        std::set<unsigned int>::const_iterator iter;
        double totalSqDist = 0.0;
        if (queryTracklet->indices.size() == 0) {
            throw LSST_EXCEPT(BadParameterException, 
                              "getAverageSqDist called with tracklet of length 0\n");
        }
        for (iter = queryTracklet->indices.begin(); iter != queryTracklet->indices.end(); iter++) {
            totalSqDist += getSqDist(slopeAndOffsetRA, slopeAndOffsetDec, &((*allDets)[*iter]), timeOffset);
        }
        return totalSqDist / queryTracklet->indices.size();
    }





    bool trackletsAreCompatible(const std::vector<MopsDetection> *detections, 
                                const Tracklet &t1, const Tracklet &t2) {
        std::set<unsigned int> uniqueIndices; 
        std::set<unsigned int>::iterator indexIter;
        std::set<unsigned int>::iterator uniqueIndicesIter;
        /* get a list of all detections in either tracklet, then get a list of
         * observation times for each tracklet, making sure there are no repeats
         * in the *observation* times.*/
        std::set_union(t1.indices.begin(), t1.indices.end(), t2.indices.begin(), t2.indices.end(),
                       std::inserter(uniqueIndices, uniqueIndices.begin()));

        std::set<double> representedMJDs;

        for (uniqueIndicesIter = uniqueIndices.begin(); uniqueIndicesIter != uniqueIndices.end(); 
             uniqueIndicesIter++) {
            double thisMJD = (*detections)[*uniqueIndicesIter].getEpochMJD();
            if (representedMJDs.insert(thisMJD).second == false) {
                /* insert failed because MJD is already there.*/
                return false;
            };
        }
        /* if we got this far without returning false, then all unique
         * detections came from separate times */
        return true;
    }



        
        void collapse(Tracklet &t1, Tracklet& t2) {
        t1.isCollapsed = true;
        t2.isCollapsed = true;
        std::set<unsigned int> uniqueIndices; 
        std::set_union(t1.indices.begin(), t1.indices.end(), t2.indices.begin(), t2.indices.end(),
                       std::inserter(uniqueIndices, uniqueIndices.begin()));
        t2.indices = uniqueIndices;
    }




    Tracklet unionTracklets(Tracklet t1, Tracklet t2) {
        Tracklet toRet;
        std::set<unsigned int> uniqueIndices; 
        std::set_union(t1.indices.begin(), t1.indices.end(), t2.indices.begin(), t2.indices.end(),
                       std::inserter(uniqueIndices, uniqueIndices.begin()));
        toRet.indices = uniqueIndices;
        toRet.isCollapsed = false;
        return toRet;
    }
    



  
    /*
     * given vector detections and pairs, a vector of vectors of indices into
     * detections, "collapse" together highly similar tracklets (where each
     * element of "pairs" describes a tracklet, a collection of detections)
     * put the resulting tracklets into collapsedPairs.
     */
    void doCollapsingPopulateOutputVector(
        const std::vector<MopsDetection> * detections, 
        std::vector<Tracklet> &pairs,
        std::vector<double> tolerances, 
        std::vector<Tracklet> &collapsedPairs,       
        bool useMinimumRMS, bool useBestFit, 
        bool useRMSFilt,
        double maxRMS, bool beVerbose) {


      time_t linkingStart = time(NULL);
        /* each t in trackletsForTree maps tracklet physical parameters (RA0,
         * Dec0, angle, vel.) to an index into pairs. */
        std::vector<PointAndValue <unsigned int> >
            trackletsForTree;

        std::vector<GeometryType> geometryTypes(4);
        /* RA0, Dec0, and angle are all degree measures along [0,360).
	   Velocity is purely euclidean.*/
        geometryTypes[0] = CIRCULAR_DEGREES;
        geometryTypes[1] = CIRCULAR_DEGREES;
        geometryTypes[2] = CIRCULAR_DEGREES;
        geometryTypes[3] = EUCLIDEAN;

        if (beVerbose) {
            std::cout << "Extrapolating linear movement functions for each tracklet." << std::endl;
        }
        populateTrackletsForTreeVector(detections, &pairs, trackletsForTree);
        if (beVerbose) {
            std::cout << "done." << std::endl;
            
            std::cout << "Building KDTree of all tracklets.." << std::endl;
        }
        KDTree<unsigned int> searchTree(trackletsForTree, 4, MAX_LEAF_SIZE);       
        if (beVerbose) {
            std::cout << "done." << std::endl;
            std::cout << "Doing many, many tree queries and collapses..." << std::endl;
        }

    int nthreads, tid;
    char* chunkSizeStr = NULL;
    int chunkSize = 4096;
    chunkSizeStr = getenv ("CHUNK_SIZE");
    if (chunkSizeStr!=NULL)
      { chunkSize = atoi(chunkSizeStr); }

#pragma omp parallel private(nthreads, tid)
    {
      nthreads = omp_get_num_threads();
      tid = omp_get_thread_num();
      if(tid == 0) {
	std::cout << "Number of threads " << nthreads << std::endl;
	std::cout << "Chunk size: " << chunkSize << std::endl;
      }
    }

#pragma omp parallel for schedule(dynamic, chunkSize)
        for (unsigned int ti = 0; ti < trackletsForTree.size(); ti++) {
            
            PointAndValue<unsigned int>* curTracklet = &(trackletsForTree[ti]);
            /* don't collapse a given tracklet twice */
            if (pairs[curTracklet->getValue()].isCollapsed == false) {
                
                /* create a new tracklet which will be output. it is marked as
                   collapsed already, and so will be the current tracklet from
                   pairs.  This way we won't bother trying to collapse this tracklet again - 
                   if we don't get it now, it won't happen later, either. */
                Tracklet newTracklet;
                collapse(pairs[curTracklet->getValue()], newTracklet);

                /* find all similar tracklets */
                std::vector<PointAndValue<unsigned int> > queryResults = 
                    searchTree.hyperRectangleSearch((*curTracklet).getPoint(), 
                                                    tolerances, 
                                                    geometryTypes);
                


                std::vector<PointAndValue<unsigned int> >::iterator trackletIter;
                std::vector<PointAndValue<unsigned int> >::iterator similarTrackletIter;
                
                if (useMinimumRMS) {
                    bool done = false;
                    /* until no work is done: check the RMS of the current
                     * tracklet combined with the each similar tracklet.  Choose
                     * the "best" option and collapse that into the current
                     * tracklet. repeat.
                     */
                    while (done == false) {
                        bool foundOne = false;
                        double bestMatchRMS = 1337;
                        unsigned int bestMatchID = 1337;
                        for (similarTrackletIter = queryResults.begin();
                             similarTrackletIter != queryResults.end();
                             similarTrackletIter++) {
                            /* try combining the current tracklet with each
                               similar tracklet and getting an RMS value. */
                            unsigned int similarTrackletID = similarTrackletIter->getValue();			    
			    if ((pairs[similarTrackletID].isCollapsed == false) && 
                                (trackletsAreCompatible(detections, pairs[similarTrackletID], newTracklet))) {
                                Tracklet tmp = unionTracklets(newTracklet, pairs[similarTrackletID]);
                                double tmpRMS = rmsForTracklet(tmp, detections);
                               if ((useRMSFilt == false) || 
                                    (tmpRMS <= maxRMS)) {
                                    if ((foundOne == false) || (bestMatchRMS > tmpRMS)) {
                                        foundOne = true;
                                        bestMatchRMS = tmpRMS;
                                        bestMatchID = similarTrackletID;
                                    }
                                }
                            }
                        }
                        if (foundOne == true) {
                            collapse(pairs[bestMatchID], newTracklet);
                        }
                        else {
                            /* found no allowable matches*/
                            done = true;
                        }
                    }
                }
                else if (useBestFit) {
                    bool done = false;
                    /* until no work is done: get the best-fit equation for the
                     * current concept of the tracklet.  Find the legal similar
                     * tracklet whose detections are closest to the current
                     * conception of the line (if one exists).  Collapse that
                     * tracklet into the current one.  Repeat.
                     */
                    while (done == false) {
                        bool foundOne = false;
                        double bestMatchAvSqDist = 1337;
                        unsigned int bestMatchID = 1337;
                        std::vector <double> currentRAFunc;
                        std::vector <double> currentDecFunc;
                        std::vector <MopsDetection> trackletDets;
                        std::set<unsigned int>::iterator detIter;
                        for (detIter = newTracklet.indices.begin(); 
                             detIter != newTracklet.indices.end(); 
                             detIter++) {
                            trackletDets.push_back((*detections)[*detIter]);
                        }
                        double t0 = (*detections)[*newTracklet.indices.begin()].getEpochMJD();
                        leastSquaresSolveForRADecLinear(&trackletDets, currentRAFunc, 
                                                                    currentDecFunc, t0);

                        for (similarTrackletIter = queryResults.begin();
                             similarTrackletIter != queryResults.end();
                             similarTrackletIter++) {
                            /* find the similar tracklet closest to the current line. */
                            unsigned int similarTrackletID = similarTrackletIter->getValue();
                            const Tracklet* similarTracklet = &pairs[similarTrackletID];
                            std::map<unsigned int, double> detIDToSqDist;
                            if ((pairs[similarTrackletID].isCollapsed == false) && 
                                (trackletsAreCompatible(detections, *similarTracklet, newTracklet))) {
                                double netSqDist = 0;
                                unsigned int newDets = 1;
                                for (std::set<unsigned int>::const_iterator sIter=similarTracklet->indices.begin();
                                     sIter != similarTracklet->indices.end(); sIter++) {
                                    if (newTracklet.indices.find(*sIter) == newTracklet.indices.end()) {
                                        // this detection is not already in our current tracklet
                                        newDets++;
                                        if (detIDToSqDist.find(*sIter) != detIDToSqDist.end()) {
                                            // look it up
                                            netSqDist += detIDToSqDist[*sIter];
                                        }
                                        else {
                                            //calculate it and save it
                                            double newSqDist = getSqDist(currentRAFunc, currentDecFunc, 
                                                                         &((*detections)[*sIter]), t0);
                                            detIDToSqDist.insert(std::make_pair(*sIter, newSqDist));
                                            netSqDist += newSqDist;
                                        }
                                    }
                                }
                                double avSqDist = netSqDist / newDets;
                                if ((foundOne == false) || (avSqDist < bestMatchAvSqDist)) {
                                    foundOne = true;
                                    bestMatchAvSqDist = avSqDist;
                                    bestMatchID = similarTrackletID;
                                }
                            }
                        }
                        if (foundOne == true) {
                            // if newTracklet + best match has higher RMS than filter allows, we're done!
                            Tracklet tmp = unionTracklets(newTracklet, pairs[bestMatchID]);
                            double tmpRMS = rmsForTracklet(tmp, detections);
                            if ((useRMSFilt == true) && (tmpRMS >  maxRMS)) {
                                done = true;
                            }
                            else { /* we got a result, and it was legal */
                                collapse(pairs[bestMatchID], newTracklet);
                            }
                        }
                        else {
                            /* found no allowable matches */
                            done = true;
                        }
                    }
                }
                else {
                    /* be greedy - try tracklets without much discriminiation */
                    for (unsigned int ri = 0; ri <  queryResults.size(); ri++) {
                        /* if tracklet is similar, has not already been collapsed,
                           not == query tracklet, and compatible with this tracklet
                           so far, then go ahead and collapse them together
                           greedily. */
                        Tracklet* similarTracklet = &(pairs[queryResults[ri].getValue()]);
                        if ((queryResults[ri].getValue() != curTracklet->getValue())
                            &&
                            (similarTracklet->isCollapsed == false)
                            && 
                            (trackletsAreCompatible(detections, *similarTracklet, newTracklet))) {
                            bool collapseIsLegal = true;
                            if (useRMSFilt) {
                                /* check that this is a 'good enough' fit to use. */
                                Tracklet tmp = newTracklet;
                                /* subtly abuse the 'collapse' function as a union operation */
                                collapse(*similarTracklet, newTracklet);
                                if (rmsForTracklet(tmp, detections) > maxRMS) {
                                    collapseIsLegal = false;
                                }
                            }
                            if (collapseIsLegal) {
                                collapse(*similarTracklet, newTracklet);
                            }
                        }
                    }                        
                }

                /* this tracklet is valid output. */
#pragma omp critical(writeOutput) 
                {
                    collapsedPairs.push_back(newTracklet);
                }                
            } /* end 'if (pairs[curTracklet->getValue().isCollapsed == false) */
        } /* end 'for (curTracklet in trackletsForTree... ) */
	std::cout << "Linking took " << timeElapsed(linkingStart) << " sec." << std::endl;
        
    }
    
















  /*
   * set define t = 0 as normalTime, so all times in the detections
   * are offsets thereof.  Then find best-fit function for mapping
   * time to RA and Dec (in deg).  Extrapolate RA, Dec at t=0; also
   * find the angle above RA = 0, Dec=t and the velocity in deg/day of
   * the tracklet.  These three will be our "point" in space.
   * 
   * ASSUME that motionVector was allocated with 4 slots.
   */
    void parameterize(const std::vector<MopsDetection> *trackletDets,
                      std::vector<double> &motionVector,
                      double normalTime) {
        std::vector<double> RASlopeAndOffset;
        std::vector<double> DecSlopeAndOffset;
        

        // physical params 0, 1 are initial RA, initial Dec at normalTime.
        leastSquaresSolveForRADecLinear(trackletDets, RASlopeAndOffset, DecSlopeAndOffset,
                                                    normalTime);
        motionVector[0] = convertToStandardDegrees(RASlopeAndOffset[1]);
        motionVector[1] = convertToStandardDegrees(DecSlopeAndOffset[1]);
        double RAv = RASlopeAndOffset[0];
        double Decv = DecSlopeAndOffset[0];
        double velocity = sqrt(pow(RAv, 2) + pow(Decv, 2));
        if (velocity != 0) {
            motionVector[3] = velocity;
            double sineOfTheta = Decv/velocity;
            double theta = asin(sineOfTheta) * (360./(2*M_PI));
            motionVector[2] = convertToStandardDegrees(theta);
        }
        else {
            motionVector[3] = 0;
            motionVector[2] = 0;
        }
    }




    /* returns average of first obs time + last obs time */
    double getMidPointTime(const std::vector<MopsDetection> *detections) {
        // these are reasonable extremes.
        double firstTime = IMPOSSIBLY_LATE_MJD;
        double lastTime = IMPOSSIBLY_EARLY_MJD;
        double midPointTime = 0.0;
        for (unsigned int i = 0; i < (*detections).size(); i++) {
            double curMJD = (*detections)[i].getEpochMJD();
            if (curMJD > lastTime) {
                lastTime = curMJD;
            }
            if (curMJD < firstTime) {
                firstTime = curMJD;
            }
        }
	midPointTime = (firstTime + lastTime) / 2.0;

        return midPointTime;
        
    }



    void populateTrackletsForTreeVector(const std::vector<MopsDetection> *detections,
                                        const std::vector<Tracklet> * tracklets,
                                        std::vector<PointAndValue <unsigned int> >
                                        &trackletsForTree) {
        
        double midPointTime = getMidPointTime(detections);
        std::vector<Tracklet>::const_iterator trackletIter;
        std::set<unsigned int>::iterator indicesIter;

        unsigned int i = 0;
        for (trackletIter = (*tracklets).begin(); trackletIter != (*tracklets).end(); trackletIter++) {
            PointAndValue <unsigned int> curTracklet;
            std::vector<double> motionVector(4); /* will hold RA0, Dec0, angle, velocity */
            std::vector<MopsDetection> trackletDets;

            for (indicesIter = trackletIter->indices.begin(); 
                 indicesIter != trackletIter->indices.end();
                 indicesIter++) {
                trackletDets.push_back((*detections)[*indicesIter]);
            }

            parameterize(&trackletDets, motionVector, midPointTime);
            /* for the KDTree, use physical params as the spatial
               'point' and make sure we get a reference to the current
               index into pairs as the associated value */
            curTracklet.setPoint(motionVector);
            /*each PointAndValue for the tree has Value which is index into
             * tracklets vector of corresponding tracklet*/
            curTracklet.setValue(i);
            trackletsForTree.push_back(curTracklet);
            i++;
        }
    }



    }} // close lsst::mops



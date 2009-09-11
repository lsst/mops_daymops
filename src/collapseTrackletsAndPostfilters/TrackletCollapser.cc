// -*- LSST-C++ -*-
/* jonathan myers */


#include <cmath>
#include <sstream>

#include <gsl/gsl_fit.h>

#include "TrackletCollapser.h"
#include "removeSubsets.h"
#include "rmsLineFit.h"

// these just need to be greater or less than any date we'll ever see...
#define IMPOSSIBLY_EARLY_MJD -1.0E16
#define IMPOSSIBLY_LATE_MJD 1.0E16

// hard to guess an ideal number for this - could use some benchmarks.
#define MAX_LEAF_SIZE 50 


namespace ctExcept = collapseTracklets::exceptions;

namespace collapseTracklets {

    // declare some helper functions - these don't need to be exposed.

    void populateTrackletsForTreeVector(const std::vector<Detection> *detections,
                                        const std::vector<Tracklet> *tracklets,
                                        std::vector<KDTree::PointAndValue <unsigned int> >
                                        &trackletsForTree);




    /* returns true iff two tracklets are 'compatible', i.e. the set of their
       combined unique detections contains no detections with the same
       observation time. NB: if t1 is detections AB, and t2 is detections AC,
       then they ARE compatible if B, C have different times (because their
       combined detections are ABC). */

    bool trackletsAreCompatible(const std::vector<Detection> * detections,
                                Tracklet t1, Tracklet t2);



    /* returns the average of the squared distances of each detection in queryTracklet from
     * the line described by slopeAndOffsetRA, slopeAndOffsetDec. */ 
    double getAverageSqDist(std::vector<double> slopeAndOffsetRA, std::vector<double> slopeAndOffsetDec,
                            const std::vector<Detection>* allDets, const Tracklet* queryTracklet, 
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



    bool TrackletCollapser::isSane(unsigned int detsSize, const std::vector<Tracklet> *pairs) {
        std::vector<Tracklet>::const_iterator tIter;
        std::set<unsigned int>::const_iterator iIter;
        for (tIter = pairs->begin(); tIter != pairs->end(); tIter++) {
            for (iIter = tIter->indices.begin(); iIter != tIter->indices.end(); 
                 iIter++) {
                if (*iIter >= detsSize) {
                    return false;
                }
            }
        }
        return true;
    }




    
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
                     const Detection* det, double timeOffset) {
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
                            const std::vector<Detection>* allDets, const Tracklet* queryTracklet, 
                            double timeOffset) {
        std::set<unsigned int>::const_iterator iter;
        double totalSqDist = 0.0;
        if (queryTracklet->indices.size() == 0) {
            throw LSST_EXCEPT(ctExcept::BadParameterException, 
                              "getAverageSqDist called with tracklet of length 0\n");
        }
        for (iter = queryTracklet->indices.begin(); iter != queryTracklet->indices.end(); iter++) {
            totalSqDist += getSqDist(slopeAndOffsetRA, slopeAndOffsetDec, &((*allDets)[*iter]), timeOffset);
        }
        return totalSqDist / queryTracklet->indices.size();
    }





    bool trackletsAreCompatible(const std::vector<Detection> *detections, Tracklet t1, Tracklet t2) {
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




    void TrackletCollapser::collapse(Tracklet &t1, Tracklet& t2) {
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
    void TrackletCollapser::doCollapsingPopulateOutputVector(
        const std::vector<Detection> * detections, 
        std::vector<Tracklet> &pairs,
        std::vector<double> tolerances, 
        std::vector<Tracklet> &collapsedPairs,       
        bool useMinimumRMS, bool useBestFit, 
        bool useRMSFilt,
        double maxRMSm, double maxRMSb, bool beVerbose) {

        /* each t in trackletsForTree maps tracklet physical parameters (RA0,
         * Dec0, angle, vel.) to an index into pairs. */
        std::vector<KDTree::PointAndValue <unsigned int> >
            trackletsForTree;
        std::vector<KDTree::PointAndValue<unsigned int> >::iterator trackletIter;
        std::vector<KDTree::PointAndValue<unsigned int> >::iterator similarTrackletIter;

        std::vector<KDTree::Common::GeometryType> geometryTypes(4);
        /* RA0, Dec0, and angle are all degree measures along [0,360).
	   Velocity is purely euclidean.*/
        geometryTypes[0] = KDTree::Common::CIRCULAR_DEGREES;
        geometryTypes[1] = KDTree::Common::CIRCULAR_DEGREES;
        geometryTypes[2] = KDTree::Common::CIRCULAR_DEGREES;
        geometryTypes[3] = KDTree::Common::EUCLIDEAN;

        if (beVerbose) {
            std::cout << "Extrapolating linear movement functions for each tracklet." << std::endl;
        }
        populateTrackletsForTreeVector(detections, &pairs, trackletsForTree);
        if (beVerbose) {
            std::cout << "done." << std::endl;
            
            std::cout << "Building KDTree of all tracklets.." << std::endl;
        }
        KDTree::KDTree<unsigned int> searchTree(trackletsForTree, 4, MAX_LEAF_SIZE);       
        if (beVerbose) {
            std::cout << "done." << std::endl;
            std::cout << "Doing many, many tree queries and collapses..." << std::endl;
        }
        unsigned int trackletCount = 0;
        for (trackletIter = trackletsForTree.begin(); 
             trackletIter != trackletsForTree.end();
             trackletIter++) {
            trackletCount++;
            /* don't collapse a given tracklet twice */
            if (pairs[trackletIter->getValue()].isCollapsed == false) {
                
                /* create a new tracklet which will be output. it is marked as
                   collapsed already, and so will be the current tracklet from
                   pairs.  This way we won't bother trying to collapse this tracklet again - 
                   if we don't get it now, it won't happen later, either. */
                Tracklet newTracklet;
                collapse(pairs[trackletIter->getValue()], newTracklet);

                /* find all similar tracklets */
                std::vector<KDTree::PointAndValue<unsigned int> > queryResults = 
                    searchTree.hyperRectangleSearch((*trackletIter).getPoint(), 
                                                    tolerances, 
                                                    geometryTypes);
                


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
                                double tmpRMS = rmsLineFit::rmsForTracklet(tmp, detections);
                                double trackletAvMag = rmsLineFit::getAverageMagnitude(tmp, detections);  
                               if ((useRMSFilt == false) || 
                                    (tmpRMS <= trackletAvMag * maxRMSm + maxRMSb)) {
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
                        std::vector <Detection> trackletDets;
                        std::set<unsigned int>::iterator detIter;
                        for (detIter = newTracklet.indices.begin(); 
                             detIter != newTracklet.indices.end(); 
                             detIter++) {
                            trackletDets.push_back((*detections)[*detIter]);
                        }
                        double t0 = (*detections)[*newTracklet.indices.begin()].getEpochMJD();
                        leastSquaresSolveForRADecLinear(&trackletDets, currentRAFunc, currentDecFunc, t0);

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
                            double tmpRMS = rmsLineFit::rmsForTracklet(tmp, detections);
                            double trackletAvMag = rmsLineFit::getAverageMagnitude(tmp, detections); 
                            if ((useRMSFilt == true) && 
                                (tmpRMS > trackletAvMag * maxRMSm + maxRMSb)) {
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
                    for (similarTrackletIter = queryResults.begin();
                         similarTrackletIter != queryResults.end();
                         similarTrackletIter++) {
                        /* if tracklet is similar, has not already been collapsed,
                           not == query tracklet, and compatible with this tracklet
                           so far, then go ahead and collapse them together
                           greedily. */
                        Tracklet similarTracklet = pairs[similarTrackletIter->getValue()];
                        if ((similarTrackletIter->getValue() != trackletIter->getValue()) 
                            &&
                            (similarTracklet.isCollapsed == false)
                            && 
                            (trackletsAreCompatible(detections, similarTracklet, newTracklet))) {
                            bool collapseIsLegal = true;
                            if (useRMSFilt) {
                                /* check that this is a 'good enough' fit to use. */
                                Tracklet tmp = newTracklet;
                                /* subtly abuse the 'collapse' function as a union operation */
                                collapse(pairs[similarTrackletIter->getValue()], newTracklet);
                                if (rmsLineFit::rmsForTracklet(tmp, detections) >
                                    rmsLineFit::getAverageMagnitude(tmp, detections) * maxRMSm + maxRMSb) {
                                    collapseIsLegal = false;
                                }
                            }
                            if (collapseIsLegal) {
                                collapse(pairs[similarTrackletIter->getValue()], newTracklet);
                            }
                        }
                    }                        
                }

                /* this tracklet is valid output. */
                collapsedPairs.push_back(newTracklet);
            } /* end 'if (pairs[trackletIter->getValue().isCollapsed == false) */
        
        } /* end 'for (trackletIter in trackletsForTree... ) */

        /* temporary sanity check */

        for (std::vector<Tracklet>::iterator tIter = pairs.begin();
             tIter != pairs.end();
             tIter++) {
            if (tIter->isCollapsed == false) {
                std::cerr << "got a tracklet which was not collapsed!" << std::endl;
            }
        }
    }
    




    void TrackletCollapser::writeTrackletsToOutFile(const std::vector<Tracklet> * tracklets, std::ofstream &outFile){
        std::vector<Tracklet>::const_iterator trackletIter;
        std::set<unsigned int>::iterator indicesIter;
        for (trackletIter = tracklets->begin(); trackletIter != tracklets->end(); trackletIter++) {
            for (indicesIter = trackletIter->indices.begin(); 
                 indicesIter != trackletIter->indices.end();
                 indicesIter++) {
                outFile << *indicesIter;
                outFile << " ";
            }
            outFile << std::endl;
        }
    }







    void make180To360Negative(std::vector<double> &v) {
        std::vector<double>::iterator iter;
        for (iter = v.begin(); iter != v.end(); iter++){
            if (*iter > 180) {
                *iter = *iter - 360;
            }
        }
    }



    void TrackletCollapser::leastSquaresSolveForRADecLinear(const std::vector <Detection> *trackletDets,
                                                            std::vector<double> &RASlopeAndOffsetOut,
                                                            std::vector<double> &DecSlopeAndOffsetOut,
                                                            double timeOffset) {
        unsigned int numDets = (*trackletDets).size();
        
        std::vector<double> MJDs(numDets);
        std::vector<double> RAs(numDets);
        std::vector<double> Decs(numDets);

        if ((RASlopeAndOffsetOut.size() != 0) || (DecSlopeAndOffsetOut.size() != 0)) {
            throw LSST_EXCEPT(ctExcept::BadParameterException, 
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
            throw LSST_EXCEPT(ctExcept::ProgrammerErrorException, 
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
            throw LSST_EXCEPT(ctExcept::GSLException, "EE: gsl_fit_linear unexpectedly returned error.\n");
        }

        free(arrayRAs);
        free(arrayDecs);
        free(arrayMJDs);

    }



  /*
   * set define t = 0 as normalTime, so all times in the detections
   * are offsets thereof.  Then find best-fit function for mapping
   * time to RA and Dec (in deg).  Extrapolate RA, Dec at t=0; also
   * find the angle above RA = 0, Dec=t and the velocity in deg/day of
   * the tracklet.  These three will be our "point" in space.
   * 
   * ASSUME that physicalParams was allocated with 4 slots.
   */
    void TrackletCollapser::setPhysicalParamsVector(const std::vector<Detection> *trackletDets,
                                                    std::vector<double> &physicalParams,
                                                    double normalTime) {
        std::vector<double> RASlopeAndOffset;
        std::vector<double> DecSlopeAndOffset;
        

        // physical params 0, 1 are initial RA, initial Dec at normalTime.
        leastSquaresSolveForRADecLinear(trackletDets, RASlopeAndOffset, DecSlopeAndOffset,
                                        normalTime);
        physicalParams[0] = KDTree::Common::convertToStandardDegrees(RASlopeAndOffset[1]);
        physicalParams[1] = KDTree::Common::convertToStandardDegrees(DecSlopeAndOffset[1]);
        double RAv = RASlopeAndOffset[0];
        double Decv = DecSlopeAndOffset[0];
        double velocity = sqrt(pow(RAv, 2) + pow(Decv, 2));
        physicalParams[3] = velocity;
        double sineOfTheta = Decv/velocity;
        double theta = asin(sineOfTheta) * (360./(2*M_PI));
        physicalParams[2] = KDTree::Common::convertToStandardDegrees(theta);
    }




    /* returns average of first obs time + last obs time */
    double getMidPointTime(const std::vector<Detection> *detections) {
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




    void TrackletCollapser::populateTrackletsForTreeVector(const std::vector<Detection> *detections,
                                                           const std::vector<Tracklet> * tracklets,
                                                           std::vector<KDTree::PointAndValue <unsigned int> >
                                                           &trackletsForTree) {

        double midPointTime = getMidPointTime(detections);
        std::vector<Tracklet>::const_iterator trackletIter;
        std::set<unsigned int>::iterator indicesIter;

        unsigned int i = 0;
        for (trackletIter = (*tracklets).begin(); trackletIter != (*tracklets).end(); trackletIter++) {
            KDTree::PointAndValue <unsigned int> curTracklet;
            std::vector<double> physicalParams(4); /* will hold RA0, Dec0, angle, velocity */
            std::vector<Detection> trackletDets;

            for (indicesIter = trackletIter->indices.begin(); 
                 indicesIter != trackletIter->indices.end();
                 indicesIter++) {
                trackletDets.push_back((*detections)[*indicesIter]);
            }

            setPhysicalParamsVector(&trackletDets, physicalParams, midPointTime);
	    /* for the KDTree, use physical params as the spatial
	       'point' and make sure we get a reference to the current
	       index into pairs as the associated value */
	    curTracklet.setPoint(physicalParams);
            /*each PointAndValue for the tree has Value which is index into
             * tracklets vector of corresponding tracklet*/
	    curTracklet.setValue(i);
	    trackletsForTree.push_back(curTracklet);
            i++;
        }
    }





    void TrackletCollapser::populateDetVectorFromFile(std::ifstream &detsFile, std::vector <Detection> &myDets) {
        std::string line;
        Detection tmpDet;
        line.clear();
        std::getline(detsFile, line);
        while (detsFile.fail() == false) {
            tmpDet.fromMITIString(line);
            myDets.push_back(tmpDet);
            line.clear();
            std::getline(detsFile, line);
        }
    }






    void TrackletCollapser::populatePairsVectorFromFile(std::ifstream &pairsFile, 
                                                        std::vector<Tracklet> &pairsVector) {

        Tracklet tmpPair;
        std::string line;
        int tmpInt = -1;
        tmpPair.indices.clear();
        tmpPair.isCollapsed = false;
        line.clear();
        std::getline(pairsFile, line);
        while (pairsFile.fail() == false) {
            std::istringstream ss(line);
            ss.exceptions(std::ifstream::failbit | std::ifstream::badbit);    
            while (!ss.eof()) {
                try {
                    ss >> tmpInt;
                    ss >> std::ws;
                }
                catch (...) {
                    throw LSST_EXCEPT(ctExcept::InputFileFormatErrorException, "Improperly-formatted pairs file.\n");
                }
                tmpPair.indices.insert(tmpInt);
            }
            if (tmpPair.indices.size() < 2) {
                throw LSST_EXCEPT(ctExcept::InputFileFormatErrorException, "EE: CollapseTracklets: pairs in pairs file must be length >= 2!\n");
            }
            pairsVector.push_back(tmpPair);
            tmpPair.indices.clear();
            tmpPair.isCollapsed = false;
            line.clear();
            std::getline(pairsFile, line);   
        }
    }
  


    /* these overloaded versions are mainly for SWIG, since Python file objects != fstreams */
    void TrackletCollapser::writeTrackletsToOutFile(const std::vector<Tracklet> * tracklets, std::string outFileName)
    {
        std::ofstream outFile;
        outFile.open(outFileName.c_str());
        if (!outFile.is_open()) {
            throw LSST_EXCEPT(ctExcept::FileException,
                              "Failed to open output file " + outFileName + " - do you have permission?\n");
        }
        writeTrackletsToOutFile(tracklets, outFile);
    }
    
    void TrackletCollapser::populateDetVectorFromFile(std::string detsFileName, std::vector <Detection> &myDets)
    {
        std::ifstream detsFile;
        detsFile.open(detsFileName.c_str());
        if (!detsFile.is_open()) {
            throw LSST_EXCEPT(ctExcept::FileException,
                              "Failed to open dets file " + detsFileName + " - does this file exist?\n");
        }
        populateDetVectorFromFile(detsFile, myDets);
        
    }
    
    void TrackletCollapser::populatePairsVectorFromFile(std::string pairsFileName,
                                     std::vector <Tracklet> &pairsVector)
    {
        std::ifstream pairsFile;
        pairsFile.open(pairsFileName.c_str());
        if (!pairsFile.is_open()) {
            throw LSST_EXCEPT(ctExcept::FileException,
                              "Failed to open pairs file " + pairsFileName + " - does this file exist?\n");
        }
        populatePairsVectorFromFile(pairsFile, pairsVector);
        
    }



}



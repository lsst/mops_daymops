// -*- LSST-C++ -*-

/*
 * jmyers 7/29/08
 * 
 * collapse-tracklets is my specialized Hough transform code.  It is used for
 * quickly finding tracklets which describe similar linear motion and joining
 * them into longer tracklets.  Currently it is intended to be run from the
 * command-line, but the interface exported here should be fairly easy to port
 * into an LSST pipeline.
 *
 */



#ifndef LSST_COLLAPSE_TRACKLETS_H
#define LSST_COLLAPSE_TRACKLETS_H

#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "../KDTree.h"
#include "../Detection.h" 
#include "../Exceptions.h"
#include "../Tracklet.h"
#include "lsst/pex/exceptions.h"


namespace collapseTracklets {

    class TrackletCollapser {
    public:
        
        /* if any index in pairs is >= detsSize, return false, else return true.
         * ALWAYS do this sanity check before anything else - particularly
         * setPhysicalParams vector or similar.
         */

        bool isSane(unsigned int detsSize, const std::vector<Tracklet> *pairs);
        
        
        /* add all det indices from t1 (which are not already in t2) into t2. 
         * set t1 and t2.isCollapsed() = True.
         */
        void collapse(Tracklet &t1, Tracklet &t2);
        
        /*
          given an array of detections, solve for RA and Dec in terms of MJD on the
          assumption that both are linear.  Output vectors will look like this when done:
          
          if RA = MJD * m + b
          RASlopeAndOffsetOut = [m, b] 
          
          (and similarly for Dec.)
          
          if timeOffset is specified as non-zero, this value will be substracted from the 
          MJDs of each detection.
        */
        void leastSquaresSolveForRADecLinear(const std::vector <Detection> *trackletDets,
                                             std::vector<double> &RASlopeAndOffsetOut,
                                             std::vector<double> &DecSlopeAndOffsetOut, 
                                             double timeOffset=0.0);
        
        void writeTrackletsToOutFile(const std::vector<Tracklet> * tracklets, std::ofstream &outFile);


        /* *pairs is modified - the Tracklets will have isCollapsed set. collapsedPairs will
         * actual output data.  I.e. if pairs contains similar tracklets [1,2]  and [2,3]  they will be
         * marked as collapsed, and [1,2,3] will be added to the collapsedPairs vector.*/
        void doCollapsingPopulateOutputVector(
            const std::vector<Detection> * detections, 
            std::vector<Tracklet> &pairs,
            std::vector<double> tolerances, 
            std::vector<Tracklet> &collapsedPairs,
            bool useMinimumRMS, bool useBestFit, 
            bool useRMSFilt, double maxRMSm, double maxRMSb, bool beVerbose);
      
        void populateDetVectorFromFile(std::ifstream &detsFile, std::vector <Detection> &myDets);

        void populatePairsVectorFromFile(std::ifstream &pairsFile,
                                         std::vector <Tracklet> &pairsVector);

        void setPhysicalParamsVector(const std::vector<Detection> *trackletDets,
                                     std::vector<double> &physicalParams,
                                     double normalTime);
            

        void populateTrackletsForTreeVector(const std::vector<Detection> *detections,
                                            const std::vector<Tracklet> * tracklets,
                                            std::vector<KDTree::PointAndValue <unsigned int> >
                                            &trackletsForTree);

        /* these overloaded versions are mainly for SWIG, since Python file objects != fstreams */
        void writeTrackletsToOutFile(const std::vector<Tracklet> * tracklets, std::string outFileName);
        void populateDetVectorFromFile(std::string detsFileName, std::vector <Detection> &myDets);
        void populatePairsVectorFromFile(std::string pairsFileName,
                                         std::vector <Tracklet> &pairsVector);
            
    };
    
}

#endif

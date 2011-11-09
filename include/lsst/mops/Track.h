// -*- LSST-C++ -*-
/* jonathan myers 
   8/10/09
*/

#ifndef LSST_TRACK_H
#define LSST_TRACK_H

#include <set>
#include <vector>
#include <Eigen/Dense>
#include "gsl/gsl_cdf.h"

#include "MopsDetection.h"
#include "Tracklet.h"

/* A track is really a set of *DETECTIONS*.  We also allow the user to track the
   *tracklets* which were used to build the track, but this is just for their
   reference.

   note that a track is this is not merely the union of the detections
   associated with the tracklets associated with the; if the track contains two
   tracklets which are in 'conflict' (contain multiple detections from the same
   image time) then only one of the possible detections is chosen.

   Note that detection are stored both as diaIds as well as indices into some
   array of detections.  The diaIds are written to file and used for provenance
   as well as for computing whether two trackSets are equal or subsets of each
   other; when running linkTracklets (or other algorithms) it is usually smart
   to access the detections using indices into an array (as we can do this in
   constant time, rather than a lookup which will take at least log time).

   Use the utility function addDetection() and it will update both the
   detectionIndices and detectionDiaIds appropriately.  

   It is up to the user to keep a Track's Detection Indices associated with the
   same Detection vector throughout its lifetime; if you have two distinct
   Detection vectors (currently this has no use case) and you create a track
   using indices from one vector, then try to refer to the other vector using
   indices from the same track you will get into trouble!

*/

namespace lsst { namespace mops {

class Track {
public:

    Track();
    
    void addDetection(unsigned int detIndex, const std::vector<MopsDetection> & allDets);

    /* Add a tracklet to the track; this will add the DETECTIONS of the tracklet
       to the current track's detection set AND add the tracklet's given tracklet index
       to componentTrackletIndices.  Used in LinkTracklets to quickly add endpoint tracklets.

       DO NOT BE MISLEAD: This is not the only way in which detection/tracklets
       are added in linkTracklets. Support detections may be added without any
       associated tracklets.
    */
    void addTracklet(unsigned int trackletIndex, 
                     const Tracklet &t, 
                     const std::vector<MopsDetection> & allDets);

    const std::set<unsigned int> getComponentDetectionIndices() const;

    const std::set<unsigned int> getComponentDetectionDiaIds() const;

    double getProbChisqRa() const { return probChisqRa; }
    double getProbChisqDec() const { return probChisqDec; }
    double getFitRange() const;

    /* until this function is called, initial position, velocity and
       acceleration for the track are NOT SET.  the USER is responsible for
       calling before using predictLocationAtTime() or getBestFitQuadratic().
     */
    void calculateBestFitQuadratic(const std::vector<MopsDetection> &allDets,
                                   const bool useFullRaFit=false, std::ostream *outFile = NULL);
    
    /* use best-fit quadratic to predict location at time mjd. will return WRONG VALUES
     if calculateBestFitQuadratic has not been called.*/
    void predictLocationAtTime(const double mjd, double &ra, double &dec) const;
    
    /* use best-fit quadratic to predict location uncertainty at time mjd. will return WRONG VALUES
     if calculateBestFitQuadratic has not been called.*/
    void predictLocationUncertaintyAtTime(const double mjd, double &raUnc, double &decUnc) const;
    
    /* you MUST call calculateBestFitQuadratic before calling this. */
    void getBestFitQuadratic(double &epoch,
                             double &ra0, double &raV, double &raAcc,
                             double &dec0, double &decV, double &decAcc) const;
    
    /*
      the tracklets which were used to build this track, if any. Currently this
      information is not used but could be useful for debugging or investigation.
    */
    std::set<unsigned int> componentTrackletIndices;
    

    Track & operator= (const Track &other);

    bool operator==(const Track &other) const {
        bool toRet = componentDetectionDiaIds == other.componentDetectionDiaIds;
        return toRet ;
    }

    bool operator!=(const Track &other) const {
        return ! (*this == other);
    }


    int getObjectId(std::vector<MopsDetection> allDets);
        

    /* the results of this comparison are probably not meaningful to a human but
     * this operator is needed for building container classes for this class */
    bool operator<(const Track &other) const {
        bool toRet = componentDetectionDiaIds < other.componentDetectionDiaIds;
        return toRet;

    }

private:
    std::set<unsigned int> componentDetectionIndices;
    std::set<unsigned int> componentDetectionDiaIds;
    Eigen::VectorXd raFunc;
    Eigen::VectorXd decFunc;
    Eigen::MatrixXd raCov;
    Eigen::MatrixXd decCov;
    double chisqRa;
    double chisqDec;
    double probChisqRa;
    double probChisqDec;
    double epoch;
    double meanTopoCorr;
};

}} // close lsst::mops namespace


#endif

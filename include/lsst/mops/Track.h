// -*- LSST-C++ -*-
/* jonathan myers 
   8/10/09
*/

#ifndef LSST_TRACK_H
#define LSST_TRACK_H

#include <set>
#include <vector>

#include "MopsDetection.h"

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

    void addDetection(unsigned int detIndex, const std::vector<MopsDetection> & allDets);

    Track();

    const std::set<unsigned int> getComponentDetectionIndices() const;

    const std::set<unsigned int> getComponentDetectionDiaIds() const;

    /* until this function is called, initial position, velocity and acceleration
       for the track are NOT SET.  
     */
    void calculateBestFitQuadratic(const std::vector<MopsDetection> &allDets);
    
    /* use best-fit quadratic to predict location at time mjd. will return WRONG VALUES
     if calculateBestFitQuadratic has not been called.*/
    void predictLocationAtTime(const double mjd, double &ra, double &dec) const;
    
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


    /* the results of this comparison are probably not meaningful to a human but
     * this operator is needed for building container classes for this class */
    bool operator<(const Track &other) const {
        bool toRet = componentDetectionDiaIds < other.componentDetectionDiaIds;
        return toRet;
        
    }

private:
    std::set<unsigned int> componentDetectionIndices;
    std::set<unsigned int> componentDetectionDiaIds;
    std::vector<double> raFunc;
    std::vector<double> decFunc;
    double epoch;
    void  bestFit1d(const std::vector<double> &X,
                    const std::vector<double> &time,
                    std::vector<double> & res);
};

}} // close lsst::mops namespace


#endif

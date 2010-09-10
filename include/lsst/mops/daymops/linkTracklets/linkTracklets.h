// -*- LSST-C++ -*-
/* jonathan myers */

#ifndef LSST_LINKTRACKLETS_H_
#define LSST_LINKTRACKLETS_H_

#include <vector>

#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/Tracklet.h"
#include "lsst/mops/Track.h"
#include "lsst/mops/TrackSet.h"


namespace lsst {
    namespace mops {

enum trackOutputMethod { RETURN_TRACKS = 0, 
                         IDS_FILE,
                         IDS_FILE_WITH_CACHE};

/*
 * linkTrackletsConfig is a simple class which holds the many configurable
 * parameters for running linkTracklets.  The default values should be 
 * fine for LSST.
 */

class linkTrackletsConfig {
public:

    linkTrackletsConfig() 
        {   maxRAAccel = .02; // consistent with 99% of MBOs
            maxDecAccel = .02;  // consistent with 99% of MBOs
            detectionLocationErrorThresh = .0005; // .3 arcseconds = 8.3e-5 degrees
            minEndpointTimeSeparation = 2; 
            minSupportToEndpointTimeSeparation = .5;
            minSupportTracklets = 1;
            quadraticFitErrorThresh = .000; // be optimistic for now...
            minDetectionsPerTrack = 6;

            // Kubica sets vtree_thresh to .0002 degrees and it's fast and accurate enough
            // If we go with a worst-case and guess that our tracklets happen at
            // most 30 min apart, this corresponds to a velocityErrorThresh of :
            // .0002 degrees * 2 / (30 min in days) = .0192 deg/day
            velocityErrorThresh = .0192;

            leafSize=16;

            restrictTrackStartTimes = false;
            latestFirstEndpointTime = -1;
            restrictTrackEndTimes = false;
            earliestLastEndpointTime = -1;


            // options for how to write output - return tracks, or write them
            // semi-continuously to disk?
            outputMethod = RETURN_TRACKS;
            outputFile = "";
            outputBufferSize = 0;

        }

    /* acceleration terms are in degrees/(day^2)
       
       objects which accelerate faster than these maximum acceleration params
       are not reported (or considered)
     */
    double maxRAAccel;
    double maxDecAccel;

    /* 
       if you will perform repeated runs, it may be wise to look for tracks
       which start or end on a limited set of nights or times.

       If you'd like to restrict linkTracklets to finding only tracks which
       start at or before time X, set restrictTrackStarTimes to true and
       latestFirstEndpointTime to X.

       if you'd like to look only for tracks which *end* after time Y (that is,
       their final endpoint tracklet begins on time Y or later) then set
       restrictTrackEndTimes to true and set earliestLastEndpointTime to Y.
     */
    bool restrictTrackStartTimes;
    double latestFirstEndpointTime;

    bool restrictTrackEndTimes;
    double earliestLastEndpointTime;

    /* detection error thresh is the upper bound on observational error for
       Detections; this has repercussions for what Detections are included in a
       returned track as well as the behavior of the searching itself (we will
       prune off regions of space based on their location, so error threshold
       will factor in here.)
     */
    double detectionLocationErrorThresh;

    
    /* trackletVelocityErrorThresh: 
     *
     * used in tree calculations to determine the maximum feasible error of a tracklet's velocity.
     * should be >= detectionLocationErrorThresh * 2.0  / max time between any two images.
     *
     * JMYERS: No, wait, shouldn't it be detectionLocationErrorThresh * 2.0 / __MIN__ time between any two images?
     */
    double trackletVelocityErrorThresh;

    /* quadratic fit error thresh: 

       this is used to determine how much error we will accept from the
       quadratic fitting process.  If a detection is within
       detectionLocationErrorThresh + quadraticFitErrorThresh of the predicted
       location, it will be considered for addition to a track
     */
    double quadraticFitErrorThresh;


    /*
      minSupportToEndpointTimeSeparation: the minimum time between endpoints and
      support points.  normally, orbit fits need several days of observations in
      order to be useful, so we will ignore support points which happen right
      after the actual detection. This also speeds up the search a little.
    */
    double minSupportToEndpointTimeSeparation;
    /*
      min time between the first and last point of a tracklet in order to be
      considered for a track
     */
    double minEndpointTimeSeparation;
    /*
      minimum number of *support* tracklets required for a track.  note that the
      track will have two "endpoint" tracklets PLUS this many support tracklets.
     */
    double minSupportTracklets;
    
    /*
      the minimum number of detections per track
     */
    unsigned int minDetectionsPerTrack;

    /*
      velocity error thresh: this is used when calculating whether two regions
      could contain the same object.  this is the maximum believable discrepancy
      between *actual* velocity and best-fit velocity, measured in deg/day.
     */
    double velocityErrorThresh;

    /* 
       max # tracklets per each leaf node of the tree. Affects performance but
       should not affect correctness.
     */
    unsigned int leafSize;


    // outputMethod, outputFile, outputBufferSize:  these 
    // define how linktracklets writes its results.

    // if outputMethod is RETURN_TRACKLETS then actually buffer all results in
    // memory, return them in-memory.

    // if outputMethod is IDS_FILE or IDS_FILE_WITH_CACHE, then write a plain-text file, with one
    // track per line, written as a series of space-delimited Detection IDs
    // which comprise the track.  

    // if IDS_FILE, write to outputBufferSize all at once, when finished.

    // if IDS_FILE_WITH_CACHE, Write to a file named by outputFile, buffering
    // outputBufferSize results between writes.
    trackOutputMethod outputMethod;
    std::string outputFile;
    unsigned int outputBufferSize;

};








/* queryTracklets are non-const because we set their velocityRA and velocityDec fields. 
   otherwise queryTracklets will not be changed. */
TrackSet* linkTracklets(const std::vector<MopsDetection> &allDetections,
                        std::vector<Tracklet> &queryTracklets,
                        linkTrackletsConfig searchConfig);








//only declared for unit testing
void getBestFitVelocityAndAcceleration(std::vector<double> positions, const std::vector<double>&times,
                                       double & velocity, double &acceleration, double &position0);

// note that position and velocity will be MODIFIED. 
void modifyWithAcceleration(double &position, double &velocity, 
                            double acceleration, double time);



    }} // close lsst::mops

#endif 

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


    // Default values set in this constructor. See each variable's declaration below for
    // some descriptions of what each parameter means. 

    linkTrackletsConfig() 
        {   
            maxRAAccel = .02; // consistent with 99% of MBOs
            maxDecAccel = .02;  // consistent with 99% of MBOs


            /* 
             * Kubica has vtrees_thresh set to .0002 by default and that's what
             * we've been using.  detectionLocationErrorThresh is the same but
             * only applies to the believable error on POSITIONS (ra, dec) of
             * the KD-Tree nodes.  However we could probably set this lower - 
             * we expect most of our detections to be within .3 arcesconds of correct.
             */
            detectionLocationErrorThresh = 0.00020000; 
            minEndpointTimeSeparation = 2; 
            minSupportToEndpointTimeSeparation = .5;
            minUniqueNights = 3;
            minDetectionsPerTrack = 6;
           

            /* Kubica's default value times two.  Kubica uses .0005 deg, but
             *  then only adds the point if dist/2.0 < pred_fit, where dist is
             *  angular distance from the best-fit predicted location for the
             *  track. */
            trackAdditionThreshold = .001; 


            /* This is the square root of value we've been
             * giving Kubica (0.00000025); he takes MEAN SQUARED not ROOT mean squared as we
             * do. */           
            trackMaxRms = .0005;


            /* Kubica sets vtree_thresh to .0002 degrees and it's fast and accurate enough
             * If we go with a worst-case and guess that our tracklets happen at
             * most 30 min apart, this corresponds to a velocityErrorThresh of :
             * .0002 degrees * 2 / (30 min in days) = .0192 deg/day
             *
             * TBD: JMYERS: we may actually want to consider raising this, since
             * our min time separation is actually more like 15 minutes now. */
            
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
       Detections; this is used to decide whether bounding boxes (KD tree nodes)
       of Tracklets are compatible (could contain a track).
     */
    double detectionLocationErrorThresh;

    
    /* trackletVelocityErrorThresh: 
     *
     * used in tree calculations to determine the maximum feasible error of a
     * tracklet's velocity.  should be >= detectionLocationErrorThresh * 2.0 /
     * min time between any two images.  This is used in searching to determine
     * whether bounding boxes (KD Tree nodes) of Tracklets are compatible (could
     * contain a track)
     */

    double trackletVelocityErrorThresh;


    /* trackAdditionThreshold:
     *
     * once we have taken a pair of endpoint tracklets and fit an initial
     * quadratic curve to those tracklets, we add support detections if they are
     * within trackAdditionThreshold of the quadratic.
     */
    double trackAdditionThreshold; 

    /* trackMaxRms:
     *
     * After choosing endpoint tracklets and support detections to build a
     * track, we will calculate the overall RMS detection location error of the
     * track.  If the RMS is below trackMaxRms, then the track is added to
     * output, otherwise it is discarded.
     */
    double trackMaxRms;

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
     * minUniqueNights: 
     *
     * Tracks are not reported unless they contain detections from at least this
     * many distinct nights.
     */
    double minUniqueNights;
    
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
                        const linkTrackletsConfig &searchConfig);








//only declared for unit testing
void getBestFitVelocityAndAcceleration(std::vector<double> positions, const std::vector<double>&times,
                                       double & velocity, double &acceleration, double &position0);

// note that position and velocity will be MODIFIED. 
void modifyWithAcceleration(double &position, double &velocity, 
                            double acceleration, double time);



    }} // close lsst::mops

#endif 

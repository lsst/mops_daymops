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


/* use settings in this class to declare the verbosity of
 * linkTracklets. Useful for debugging. */
class linkTrackletsVerbositySettings 
{
public: 
    linkTrackletsVerbositySettings() {
        /* if set to true, will announce what phase of processing is
         * going on and which imag e endpoint pair is in use for
         * linking.
         */
        printStatus = false;
        printVisitCounts = false;
        printTimesByCategory = false;
        printBoundsInfo = false;
    }
    bool printStatus;
    bool printVisitCounts;
    bool printTimesByCategory;
    bool printBoundsInfo;
};




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
             * Kubica has vtrees_thresh set to .0002 by default and
             * that's what we've been using.
             * detectionLocationErrorThresh is the same but only
             * applies to the believable error on POSITIONS (ra, dec)
             * of the KD-Tree nodes.  However we could probably set
             * this lower - we expect most of our detections to be
             * within .3 arcesconds of correct.
             */
            detectionLocationErrorThresh = 0.0002; 
            // kubica sets TBT_MIN_TIME to 1. Try to replicate.
            minEndpointTimeSeparation = 1;
            minSupportToEndpointTimeSeparation = .5;
            minUniqueNights = 3;
            minDetectionsPerTrack = 6;
           

            /* Kubica uses .0005, but that's in RADIANS! this is the
             * degree equivalent. */
            trackAdditionThreshold = .028648;


            /* This is the square root of value we've been giving
             * Kubica (0.00000025); he takes MEAN SQUARED not ROOT
             * mean squared as we do. */           
            trackMaxRms = .0005;


            // Now with "sparse" KD-Trees it appears that leaf node
            // size 1 is best (see my spreadsheet on Google docs -
            // jmyers)
            leafSize=1;

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
       
       objects which accelerate faster than these maximum acceleration
       params are not reported (or considered)
     */
    double maxRAAccel;
    double maxDecAccel;


    /* 
       if you will perform repeated runs, it may be wise to look for
       tracks which start or end on a limited set of nights or times.

       If you'd like to restrict linkTracklets to finding only tracks
       which start at or before time X, set restrictTrackStarTimes to
       true and latestFirstEndpointTime to X.

       if you'd like to look only for tracks which *end* after time Y
       (that is, their final endpoint tracklet begins on time Y or
       later) then set restrictTrackEndTimes to true and set
       earliestLastEndpointTime to Y.
     */
    bool restrictTrackStartTimes;
    double latestFirstEndpointTime;

    bool restrictTrackEndTimes;
    double earliestLastEndpointTime;


    /* detection error thresh is the upper bound on observational
       error for Detections; this is used to decide whether bounding
       boxes (KD tree nodes) of Tracklets are compatible (could
       contain a track).  It is also used to derive velocity error for
       each tracklet, which trickles up to the error bounds on tree
       nodes (bounding boxes).
     */
    double detectionLocationErrorThresh;

    


    /* trackAdditionThreshold:
     *
     * once we have taken a pair of endpoint tracklets and fit an
     * initial quadratic curve to those tracklets, we add support
     * detections if they are within trackAdditionThreshold of the
     * quadratic.
     */
    double trackAdditionThreshold; 

    /* trackMaxRms:
     *
     * After choosing endpoint tracklets and support detections to
     * build a track, we will calculate the overall RMS detection
     * location error of the track.  If the RMS is below trackMaxRms,
     * then the track is added to output, otherwise it is discarded.
     */
    double trackMaxRms;

    /*
      minSupportToEndpointTimeSeparation: the minimum time between
      endpoints and support points.  normally, orbit fits need several
      days of observations in order to be useful, so we will ignore
      support points which happen right after the actual
      detection. This also speeds up the search a little.
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
     * Tracks are not reported unless they contain detections from at
     * least this many distinct nights.
     */
    double minUniqueNights;
    
    /*
      the minimum number of detections per track
     */
    unsigned int minDetectionsPerTrack;

    /* 
       max # tracklets per each leaf node of the tree. Affects
       performance but should not affect correctness.
     */
    unsigned int leafSize;


    // outputMethod, outputFile, outputBufferSize: these define how
    // linktracklets writes its results.

    // if outputMethod is RETURN_TRACKLETS then actually buffer all
    // results in memory, return them in-memory.

    // if outputMethod is IDS_FILE or IDS_FILE_WITH_CACHE, then write
    // a plain-text file, with one track per line, written as a series
    // of space-delimited Detection IDs which comprise the track.

    // if IDS_FILE, write to outputBufferSize all at once, when finished.

    // if IDS_FILE_WITH_CACHE, Write to a file named by outputFile,
    // buffering outputBufferSize results between writes.
    trackOutputMethod outputMethod;
    std::string outputFile;
    unsigned int outputBufferSize;

    linkTrackletsVerbositySettings myVerbosity;

};







/* queryTracklets are non-const because we set their velocityRA and
   velocityDec fields.  otherwise queryTracklets will not be
   changed. Detections will be recentered, though.*/
TrackSet* linkTracklets(std::vector<MopsDetection> &allDetections,
                        std::vector<Tracklet> &queryTracklets,
                        const linkTrackletsConfig &searchConfig);









// note that position and velocity will be MODIFIED. 
void modifyWithAcceleration(double &position, double &velocity, 
                            double acceleration, double time);



    }} // close lsst::mops

#endif 

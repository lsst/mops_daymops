// -*- LSST-C++ -*-
/******************************************************
 *
 * Author: Matthew Cleveland 
 * File: findTracklets.h
 *
 * heavily modified by jmyers in april of 2010, replacing
 * old args with a new findTrackletsConfig class.
 ******************************************************/

#ifndef __FINDTRACKLETS_H__
#define __FINDTRACKLETS_H__

#include <vector>

#include "lsst/mops/TrackletVector.h"
#include "lsst/mops/MopsDetection.h"


namespace lsst {
    namespace mops {



/******************************************************************************
* container class to hold arguments to findTracklets. Default constructor will
* take the liberty of setting default args for you so you only have to modify
* the fields if you're doing something special.
*******************************************************************************/

enum trackletOutputMethod { RETURN_TRACKLETS = 0, 
                            IDS_FILE,
                            IDS_FILE_WITH_CACHE};
        
class findTrackletsConfig {
public:
    
    findTrackletsConfig() 
        {   
            maxDt = .0625; // 90 minutes
            minDt = .0;
            maxV = 2.0;
            minV = 0.0;
            outputMethod = RETURN_TRACKLETS;
            outputFile = "";
            outputBufferSize = 0;
        }

    // units for these two are in days.
    // max time between two detections for any attempt to build a tracklet between them
    double maxDt;
    // min time between two detections for any attempt to build a tracklet between them
    double minDt;

    // units for these two are in deg/day.
    // maxV: maximum velocity of tracklet for which we search and return.
    double maxV;
    // minV: minimum velocity of tracklet for which we search and return.
    double minV;


    // outputMethod, outputFile, outputBufferSize:  these 
    // define how findtracklets writes its results.

    // if outputMethod is RETURN_TRACKLETS then actually buffer all results in
    // memory, return them in-memory.

    // if outputMethod is IDS_FILE or IDS_FILE_WITH_CACHE, then write a plain-text file, with one
    // tracklet per line, written as a series of space-delimited Detection IDs
    // which comprise the tracklet.  

    // if IDS_FILE, write to outputBufferSize all at once, when finished.

    // if IDS_FILE_WITH_CACHE, Write to a file named by outputFile, buffering
    // outputBufferSize results between writes.

    trackletOutputMethod outputMethod;
    std::string outputFile;
    unsigned int outputBufferSize;
};
        



/*****************************************************************
 * Main function
 * 
 * if config specifies output method of RETURN_TRACKLETS, then a pointer to a
 * TrackletVector will be returned. If other methods are specified, we simply
 * return NULL, as output is written to file.
 *****************************************************************/



TrackletVector *
findTracklets(const std::vector<MopsDetection> &allDetections, 
	      findTrackletsConfig config);

    }} // close lsst::mops

#endif

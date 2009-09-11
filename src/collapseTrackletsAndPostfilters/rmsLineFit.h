// -*- LSST-C++ -*-


/* jmyers 8/18/08
 *
 *This code is for getting the best-fit line for a series of observations
 * (i.e. a tracklet) and determining the RMS of the error for the detections.
 *
 * This also includes code for filtering out groups of tracklets which do not
 * fit the criteria.
 */


#ifndef LSST_RMSLINEFIT_H
#define LSST_RMSLINEFIT_H

#include <vector>
#include <map>

#include "../Tracklet.h"
#include "../Detection.h"
#include "../Exceptions.h"
#include "TrackletCollapser.h"


namespace rmsLineFit {

    /*
     * return the root mean squared distance of the tracklet from the line.
     * (i.e. the sqrt of the sum of the squared distances of each detection from
     * the best-fit line, all divided by the size of the tracklet).
     *
     * if perDetSqDist != NULL (the default value), then we ASSUME that
     * perDetSqDist is a pointer to an allocated vector, and we will populate it
     * s.t. mean Sq. distance (trackletDets[i],line) = perDetSqDist[i] 
     */
    double rmsForTracklet(Tracklet t, const std::vector<Detection> *detections, 
                          std::vector<double>*perDetSqDist=NULL);
    
    double getAverageMagnitude(const Tracklet t, const::std::vector<Detection>* detections);
    
/*
 * given a vector of Tracklets and the corresponding vector of Detections,
 * add to output only those tracklets for which rms < maxRMSm * av. magnitude + maxRMSm
 */
    
    void filterByLineFitAddToOutputVector(const std::vector<Tracklet> *tracklets, 
                                          const std::vector<Detection> * allDets,
                                          double maxRMSm, double maxRMSb,
                                          std::vector<Tracklet> &output);
    
    

    class TrackletPurifier {
    public:
        Tracklet purifyTracklet(const Tracklet *t, const std::vector<Detection> *allDets, 
                                double maxRMSm, double maxRMSb);

        void purifyTracklets(const std::vector<Tracklet> *trackletsVector,
                             const std::vector<Detection> *detsVector,
                             double maxRMSm, double maxRMSb, unsigned int minObs,
                             std::vector<Tracklet> &output);
    };


}

#endif

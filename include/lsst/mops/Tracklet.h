// -*- LSST-C++ -*-


/* jmyers 8/18/08
 */



#ifndef LSST_TRACKLET_H
#define LSST_TRACKLET_H
#include <set>
#include <vector>
#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/Exceptions.h"

namespace lsst {
namespace mops {

class Tracklet {
public: 
    Tracklet();
    Tracklet(std::set <unsigned int> startIndices);

    /* a tracklet is just a collection of detections.  To conserve
       memory, we don't hold copies of those detections here, we hold
       references to them.  In all my code, I represent a detection
       with its index from the Detections vector I'm working with.
       Note that these indices ARE NOT NECESSARILY the same as the
       ID's of the Detections themselves.
       
       It is TOTALLY UP TO THE USER of this class to make sure that
       they associate Tracklets with the correct Detection vector.
     */

    std::set<unsigned int> indices;
    bool isCollapsed;    
    
    void setBestFitFunctionRa(std::vector<double>);
    const std::vector<double>* getBestFitFunctionRa();
    void setBestFitFunctionDec(std::vector<double>);
    const std::vector<double>* getBestFitFunctionDec();
    
    /* return start/end times and first/last detections. does NOT
     * assume indices are assigned chronologically. DOES ASSUME,
     * however, that you send in the same vector of detections used to
     * create this tracklet! */
    double getStartTime(const std::vector<MopsDetection> &dets) const;
    MopsDetection getFirstDetection(const std::vector<MopsDetection> &dets)
        const;
    double getDeltaTime(const std::vector<MopsDetection> &dets) const;
    MopsDetection getLastDetection(const std::vector<MopsDetection> &dets)
        const;

    /*
     * The following operators consider only the indices set.  Note
     * that if you compare two tracklets with indices into different
     * detection vectors, the results are MEANINGLESS! It is up to the
     * user to ensure that when comparing tracklets t1 and t2 that if
     * t1.indices contains id X, then X in t2.indices has the same
     * meaning.
     */

    void setId(unsigned int nId) { myId = nId;}
    unsigned int getId() const { return myId; }


    Tracklet & operator= (const Tracklet &other);

    bool operator==(const Tracklet &other) const; 

    bool operator!=(const Tracklet &other) const;

    bool operator<(const Tracklet &other) const; 

private:
    std::vector<double> bestFitFunctionRa;
    std::vector<double> bestFitFunctionDec;
    void copy(const Tracklet &other);
    unsigned int myId; 

};


} } //close namespace lsst::mops

#endif

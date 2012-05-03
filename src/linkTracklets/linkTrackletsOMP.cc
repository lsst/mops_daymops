// -*- LSST-C++ -*-
/* jonathan myers */

// time headers needed for benchmarking performance
#include <ctime>
#include <iomanip>
#include <map>
#include <time.h>
#include <algorithm>
#include <omp.h>

#include "lsst/mops/rmsLineFit.h"
#include "lsst/mops/daymops/linkTracklets/linkTracklets.h"
#include "lsst/mops/Exceptions.h"
#include "lsst/mops/KDTree.h"
#include "lsst/mops/daymops/linkTracklets/TrackletTree.h"

#undef DEBUG

/* taking a queue from Kubica, it's only once per ITERATIONS_PER_SPLIT
 * calls to doLinkingRecurse that we actually split the (non-leaf)
 * support nodes. The idea is to avoid redundantly calculating whether
 * a set of support nodes is compatible.
 */
#define ITERATIONS_PER_SPLIT 0


#define POINT_RA           0
#define POINT_DEC          1
#define POINT_RA_VELOCITY   2
#define POINT_DEC_VELOCITY  3

#define uint unsigned int 

namespace lsst {
   namespace mops {








/* **********************************************************************

* A FEW SIMPLE, LOCAL CLASSES: mostly for readability and to allow the
  cache of KDTree Node, Image time, second image Time -> bounding box

* TBD: linkTracklets now uses its own special TrackletTree class so
* probably just move this stuff there?
  
*********************************************************************** */


/*
 ImageTime is a class for improving readability as well as adding
 functionality. From now on we will be giving unique IDs to each image
 time.

 This is important since we now support the use of a cache which takes
 a tree node ID and a time as a key, and you shouldn't use a
 floating-point number as a key.
 */


class ImageTime {
public:
    ImageTime() {
        MJD = -1; imgId = 0;
    }
    ImageTime(double newMJD, uint newImageId) {
        MJD = newMJD;
        imgId = newImageId;
    }
    ImageTime(const ImageTime &other) {
        MJD = other.getMJD();
        imgId = other.getImageId();
    }
    const double getMJD() const {
        return MJD;
    }
    const uint getImageId() const {
        return imgId;
    }
    void setMJD(double m) {
        MJD = m;
    }
    void setImageId(uint i) {
        imgId = i;
    }
    ImageTime & operator=(const ImageTime &rhs) {
        MJD = rhs.getMJD();
        imgId = rhs.getImageId();
        return *this;
    }
    // NB: since IDs are assigned in order of image time, would it be
    // faster if we did < based on imgId (int rather than double?)
    bool operator< (const ImageTime &other) const {
        return MJD < other.getMJD();
    }

private:
    double MJD;
    uint imgId;
};





class TreeNodeAndTime {
public:
    TreeNodeAndTime(TrackletTreeNode * tree, ImageTime i) {
        myTree = tree;
        myTime = i;
    }
    TreeNodeAndTime() {};
    TreeNodeAndTime(const TreeNodeAndTime &other) {
        myTree = other.myTree;
        myTime = other.myTime;
    }
    TreeNodeAndTime &operator=(const TreeNodeAndTime &other) {
        myTree = other.myTree;
        myTime = other.myTime;
        return *this;
    }
    TrackletTreeNode * myTree;
    ImageTime myTime;

};







/* *****************************************************
* These functions are for debugging/ diagnostics only.  
********************************************************/


std::set<uint> setUnion(const std::set<uint> &s1, 
                        const std::set<uint> &s2) 
{
    std::set<uint> toRet;
    std::set<uint>::const_iterator iter;
    for (iter = s1.begin(); iter != s1.end(); iter++) {
        toRet.insert(*iter);
    }
    for (iter = s2.begin(); iter != s2.end(); iter++) {
        toRet.insert(*iter);
    }
    return toRet;
}


std::set<uint> allDetsInTreeNode(TrackletTreeNode &t,
                                 const std::vector<MopsDetection>&allDets,
                                 const std::vector<Tracklet>&allTracklets) 
{
    std::set<uint> toRet;
    std::vector<PointAndValue<uint> >::const_iterator tIter;

    if(! t.isLeaf() ) {
        if (t.hasLeftChild()) {
            std::set<uint> childDets = allDetsInTreeNode(*t.getLeftChild(), allDets, allTracklets);
            toRet = setUnion(toRet, childDets);
        }
        if (t.hasRightChild()) {
            std::set<uint> childDets = allDetsInTreeNode(*t.getRightChild(), allDets, allTracklets);
            toRet = setUnion(toRet, childDets);
        }
    }
    else 
    {
        // t.getMyData() holds points and values. the "value" part is
        // index into allTracklets.
        for (tIter = t.getMyData()->begin();
             tIter != t.getMyData()->end();
             tIter++) {
            
            std::set<uint>::const_iterator detIter;
            for (detIter = allTracklets.at(tIter->getValue()).indices.begin();
                 detIter != allTracklets.at(tIter->getValue()).indices.end();
                 detIter++) {
                
                toRet.insert(allDets.at(*detIter).getID());
                
            }
        }
    }

    return toRet;

}





void printSet(const std::set<uint> s, std::string delimiter) 
{
    std::set<uint>::const_iterator setIter;
    for (setIter = s.begin();
         setIter != s.end();
         setIter++) {
        std::cout << *setIter << delimiter;
    }

}








void debugPrint(const TreeNodeAndTime &firstEndpoint, 
                const TreeNodeAndTime &secondEndpoint, 
                std::vector<TreeNodeAndTime> &supportNodes, 
                const std::vector<MopsDetection> &allDetections,
                const std::vector<Tracklet> &allTracklets) 
{
    std::set<uint> leftEndpointDetIds = allDetsInTreeNode(*(firstEndpoint.myTree), 
                                                          allDetections, allTracklets);
    std::set<uint> rightEndpointDetIds = allDetsInTreeNode(*(secondEndpoint.myTree),
                                                           allDetections, allTracklets);

    std::cout << "in doLinkingRecurse,       first endpoint contains     " ;
    printSet(leftEndpointDetIds, " ");
    std::cout << '\n';
    std::cout << "                            second endpoint contains    " ;
    printSet(rightEndpointDetIds, " ");
    std::cout << '\n';

    for(uint i = 0; i < supportNodes.size(); i++) {
        std::cout << " support node " << i << " contains                      ";
        std::set <uint> supDetIds = allDetsInTreeNode(*(supportNodes.at(i).myTree), 
                                                      allDetections, allTracklets);
        printSet(supDetIds, " ");
        std::cout << '\n';
    }
    

    if ((
        ((leftEndpointDetIds.find(6) != leftEndpointDetIds.end()) && 
        (rightEndpointDetIds.find(10) != rightEndpointDetIds.end())) 
        ||
        (((leftEndpointDetIds.find(12) != leftEndpointDetIds.end()) && 
          (rightEndpointDetIds.find(16) != rightEndpointDetIds.end())))
        ||
        (((leftEndpointDetIds.find(24) != leftEndpointDetIds.end()) && 
          (rightEndpointDetIds.find(28) != rightEndpointDetIds.end())))
        ||
        (((leftEndpointDetIds.find(30) != leftEndpointDetIds.end()) && 
          (rightEndpointDetIds.find(34) != rightEndpointDetIds.end())))
            ))
    {
        
        std::cout << "you'd expect to find a result here \n";
    }
    std::cout << '\n';
}




void showNumVisits(TrackletTreeNode *tree, uint &totalNodes, uint &totalVisits) 
{
    totalNodes++;
    totalVisits += tree->getNumVisits();
    std::cout << "Node " << tree->getId() << " had " << tree->getNumVisits() << " visits.\n";
    if (tree->hasLeftChild()) {
        showNumVisits(tree->getLeftChild(), totalNodes, totalVisits);
    }
    if (tree->hasRightChild()) {
        showNumVisits(tree->getRightChild(), totalNodes, totalVisits);
    }
}

















/* **********************************************************************
* LINKTRACKLETS: The actual algorithm implementation *
* **********************************************************************/



/* the final parameter is modified; it will hold Detections associated
   with the tracklet t. */

void getAllDetectionsForTracklet(
    const std::vector<MopsDetection> & allDetections,
    const Tracklet &t,
    std::vector<MopsDetection> &detectionsForTracklet) 
{
    detectionsForTracklet.clear();
    std::set<uint>::const_iterator trackletIndexIter;

    for (trackletIndexIter = t.indices.begin();
         trackletIndexIter != t.indices.end();
         trackletIndexIter++) {
        detectionsForTracklet.push_back(allDetections.at(
                                            *trackletIndexIter));
    }
}





void recenterDetections(std::vector<MopsDetection> &allDetections, 
                        const linkTrackletsConfig &searchConfig)
{
    if (allDetections.size() < 1) return;

    double centerRa = searchConfig.skyCenterRa;
    double centerDec = searchConfig.skyCenterDec;
    double minRa = centerRa;
    double maxRa = centerRa;
    double minDec = centerDec;
    double maxDec = centerDec;
    for (unsigned int i = 1; i < allDetections.size(); i++) {
        double thisRa, thisDec;
        thisRa  = allDetections[i].getRA();
        thisDec = allDetections[i].getDec();
        while (centerRa - thisRa > 180.) {
            thisRa += 360.;
        }
        while (centerRa - thisRa < -180.) {
            thisRa -= 360.;
        }
        while (centerDec - thisDec > 180.) {
            thisDec += 360.;
        }
        while (centerDec - thisDec < -180.) {
            thisDec -= 360.;
        }
        allDetections[i].setRA(thisRa);
        allDetections[i].setDec(thisDec);
        
        minRa = minOfTwo(thisRa, minRa);
        maxRa = maxOfTwo(thisRa, maxRa);

        minDec = minOfTwo(thisDec, minDec);
        maxDec = maxOfTwo(thisDec, maxDec);
    }

    if (searchConfig.myVerbosity.printBoundsInfo) {
        std::cout << "   Bounds were (in deg) R=[" << 
            minRa << ",  " << maxRa << "], D=[ " <<
            minDec << ",  " << maxDec << "]" << std::endl;
    }
    
    
}







void setTrackletVelocities(
    const std::vector<MopsDetection> &allDetections,
    std::vector<Tracklet> &queryTracklets)
    
{
    for (uint i = 0; i < queryTracklets.size(); i++) {
        Tracklet *curTracklet = &queryTracklets.at(i);
        std::vector <MopsDetection> trackletDets;
        getAllDetectionsForTracklet(allDetections, 
                                    *curTracklet, 
                                    trackletDets);

        std::vector<double> RASlopeAndOffset;
        std::vector<double> DecSlopeAndOffset;
        leastSquaresSolveForRADecLinear(&trackletDets,
                                        RASlopeAndOffset,
                                        DecSlopeAndOffset,
                                        curTracklet->getStartTime(allDetections));
        // slope and offset is reverse of p0, vel form...
        std::vector<double> raP0Vel(2);
        raP0Vel[0] = RASlopeAndOffset[1];
        raP0Vel[1] = RASlopeAndOffset[0];
        std::vector<double> decP0Vel(2);
        decP0Vel[0] = DecSlopeAndOffset[1];
        decP0Vel[1] = DecSlopeAndOffset[0];
        
        curTracklet->setBestFitFunctionRa(raP0Vel);
        curTracklet->setBestFitFunctionDec(decP0Vel);
    }

}







void makeTrackletTimeToTreeMap(
    const std::vector<MopsDetection> &allDetections,
    std::vector<Tracklet> &queryTracklets,
    std::map<ImageTime, TrackletTree > &newMap,
    const linkTrackletsConfig &myConf)
{
    bool printDebug = false;
    if (printDebug) {
        std::cout << "Sorting tracklets by \"root\" time. We have  " 
                  << allDetections.size() 
                  << " detections and " << queryTracklets.size() 
                  << " tracklets." << std::endl;

    }

    newMap.clear();

    //sort all tracklets by their first image time and assign them IDs.
    // IDs are their indices into the queryTracklets[] vector.
    
    // we use a std::vector of tracklets rather than trackletvector
    // because TrackletVectors are for large output sets, not for
    // small dynamic stuff like this. Also they have no copy
    // contructor, which is needed for map, because C++ STL doesn't
    // let you copy streams.
    std::map<double, std::vector<Tracklet> > allTrackletsByTime;

    allTrackletsByTime.clear();

    for (uint i = 0; i < queryTracklets.size(); i++) {

        double firstDetectionTime = 
            queryTracklets.at(i).getStartTime(allDetections);

        queryTracklets.at(i).setId(i);
        allTrackletsByTime[firstDetectionTime].push_back(
            queryTracklets.at(i));
    }

    if (printDebug) 
        std::cout << " got" << allTrackletsByTime.size()
                  << " image times " << std::endl;
    
    // iterate over each time/trackletVec pair and build a
    // corresponding time/KDTree pair.  note that we iterate over a
    // Map which uses MJD as key; Maps sort their data by their key,
    // so we are iterating over all image times in order.
    uint curImageId = 0;

    std::map<double, std::vector<Tracklet> >::const_iterator 
        timesIter;
    
    for (timesIter = allTrackletsByTime.begin(); 
         timesIter != allTrackletsByTime.end(); 
         timesIter++) {


        TrackletTree curTree(allDetections, 
                             timesIter->second,
                             myConf.detectionLocationErrorThresh,
                             myConf.detectionLocationErrorThresh,
                             myConf.leafSize);

        newMap[ImageTime(timesIter->first, curImageId)] = curTree;

        if (printDebug) {
            std::cout << " image time " << timesIter->first 
                      << " (with Id = " << curImageId << ") had " 
                      << timesIter->second.size() 
                      << " tracklets, generating a tree of size " 
                      << curTree.size() << std::endl;
            std::cout << "    with leaf node size = " 
                      << myConf.leafSize 
                      << ", and an average leaf size of about " 
                      << timesIter->second.size() * 2. / curTree.size() 
                      << std::endl;
        }

        curImageId++;
    }

}































/*
 * Takes a track which is ASSUMED TO HOLD exactly two tracklets; the
 * endpoint tracklets used to get a track started.  Returns true iff:
 *
 * - the first and last detection are at least minTimeSeparation apart 
 * - the best-fit accelerations are within min/max bounds
 * 
 */
// TBD: This function is a little dumb in that it looks at the time
// separation between the first and last detection, not the start time
// of the two tracklets.
bool endpointTrackletsAreCompatible(
    const std::vector<MopsDetection> & allDetections, 
    const Track &newTrack,
    const linkTrackletsConfig &searchConfig)
{
    bool allOK = true;

    double epoch, ra0, raV, raAcc;
    double dec0, decV, decAcc;
    newTrack.getBestFitQuadratic(epoch, 
                                 ra0, raV, raAcc, 
                                 dec0, decV, decAcc);
#ifdef DEBUG
    std::cout << "compatible? " << epoch << ' ' << ra0 << ' ' << raV << ' ' << raAcc
              << ' ' <<  dec0 << ' ' << decV << ' ' << decAcc  << '\n';
#endif

    if (raAcc > searchConfig.maxRAAccel) {
        allOK = false;
    }
    if (decAcc > searchConfig.maxDecAccel) {
        allOK = false;
    }

        
    if (allOK == true) {
        //check that time separation is good
        double minMJD, maxMJD;
        std::set<uint> trackDets = newTrack.getComponentDetectionIndices();
        std::set<uint>::const_iterator detIter;
        detIter = trackDets.begin();

        minMJD = allDetections.at(*detIter).getEpochMJD();
        maxMJD = minMJD;
        for (detIter = trackDets.begin();
             detIter != trackDets.end();
             detIter++) {
            double thisMJD = allDetections.at(*detIter).getEpochMJD();
            if (thisMJD < minMJD) { 
                minMJD = thisMJD; }
            if (thisMJD > maxMJD) {
                maxMJD = thisMJD;
            }
        }
        if (maxMJD - minMJD < searchConfig.minEndpointTimeSeparation) {
            allOK = false;
        }
    }

    return allOK;
}




template <class T>
bool setContains(std::set<T> s, T foo) 
{
    return (s.find(foo) != s.end());
} 



// this helper class is used just to make
// addDetectionsCloseToPredictedPositions a little more readable.

class CandidateDetection {
public:
    CandidateDetection() {
        distance = 0; 
        detId = 0;
        parentTrackletId = 0;
    }
    CandidateDetection(double nDistance, unsigned int nDetId, unsigned int nParentTrackletId) {
        distance = nDistance;
        detId = nDetId;
        parentTrackletId = nParentTrackletId;
    }
    double distance;
    unsigned int detId;
    unsigned int parentTrackletId;
};

/*
 * ASSUMES newTrack is prepared to call getBestFitQuadratic- this means
 * calculateBestFitQuadratic MUST have been called already.
 *
 * uses Track::predictLocationAtTime to find things within
 * searchConfig.trackAdditionThreshold of the predicted location.
 */
void addDetectionsCloseToPredictedPositions(
    const std::vector<MopsDetection> &allDetections, 
    const std::vector<Tracklet> &allTracklets, 
    const std::vector<uint> &candidateTrackletIds,
    Track &newTrack, 
    const linkTrackletsConfig &searchConfig)
{
    /* the scoreToIDsMap will hold detection MJDs as the key, map from
     * detection time to best-fitting candidate detection at that time
     */
    std::map<double, CandidateDetection > timeToCandidateMap;
    

    // find the best compatible detection at each unique image time
    std::vector<unsigned int>::const_iterator trackletIDIter;
    std::set<unsigned int>::const_iterator detectionIDIter;
    for (trackletIDIter = candidateTrackletIds.begin();
         trackletIDIter != candidateTrackletIds.end();
         trackletIDIter++) {
        const Tracklet * curTracklet = &allTracklets.at(*trackletIDIter);
        for (detectionIDIter =  curTracklet->indices.begin();
             detectionIDIter != curTracklet->indices.end();
             detectionIDIter++) {
            double detMjd = allDetections.at(*detectionIDIter).getEpochMJD();
            double detRa  = allDetections.at(*detectionIDIter).getRA();
            double detDec = allDetections.at(*detectionIDIter).getDec();

            double predRa, predDec;
            newTrack.predictLocationAtTime(detMjd, predRa, predDec);

            double decDistance = fabs(detDec - predDec);
#ifdef DEBUG
            std::cout << "dec dist:" << decDistance << " thresh: " << searchConfig.trackAdditionThreshold << '\n';
#endif
            if (decDistance < searchConfig.trackAdditionThreshold) {

                double distance = angularDistanceRADec_deg(detRa, detDec, predRa, predDec);
#ifdef DEBUG
                std::cout << "ang dist:" << distance << " thresh: " << searchConfig.trackAdditionThreshold << '\n';
#endif
                
                // if the detection is compatible, consider whether
                // it's the best at the image time
                if (distance < searchConfig.trackAdditionThreshold) {
                    std::map<double, CandidateDetection>::iterator candidateAtTime;
                    candidateAtTime = timeToCandidateMap.find(detMjd);
                    
                    // more crazy C++-talk for "if we have no
                    // candidate, or if this is a better candidate
                    // than the one we have, then add this as the new
                    // candidate at that time"
                    if ((candidateAtTime == timeToCandidateMap.end()) 
                        || (candidateAtTime->second.distance > distance)) {
                        CandidateDetection newCandidate(distance, 
                                                        *detectionIDIter, 
                                                        *trackletIDIter);
                        timeToCandidateMap[detMjd] = newCandidate;
                    }
                }
            }
        }
    }

    /* initialize a list of image times present in the track already. */
    std::set<double> trackMJDs;
    std::set<uint> trackDetIndicesSet = newTrack.getComponentDetectionIndices();
    std::set<unsigned int>::const_iterator trackDetectionIndices;
    for (trackDetectionIndices =  trackDetIndicesSet.begin();
         trackDetectionIndices != trackDetIndicesSet.end();
         trackDetectionIndices++) {
        trackMJDs.insert(allDetections.at(*trackDetectionIndices).getEpochMJD());        
    }
    
    /* add detections (and their parent tracklets) in order of 'score'
     * (distance from best-fit line) without adding any detections
     * from already-represented image times
     */
    std::map<double, CandidateDetection>::iterator candidatesIter;

    for (candidatesIter = timeToCandidateMap.begin(); 
         candidatesIter != timeToCandidateMap.end(); 
         candidatesIter++) {

        if (trackMJDs.find(candidatesIter->first) == trackMJDs.end()) {
            /* add this detection and tracklet */
            
            trackMJDs.insert(candidatesIter->first);
            newTrack.addDetection(candidatesIter->second.detId, 
                                  allDetections);
#ifdef DEBUG
            std::cout << "inserted detection: " << candidatesIter->second.detId << '\n';
#endif
            newTrack.componentTrackletIndices.insert(
                candidatesIter->second.parentTrackletId);
        }
    }
}




/* 
 * Return true iff the track has enough support points, and they fall
 * on enough unique nights.
 *
 * WARNING HACKISHNESS IN ACTION: see comments on how we count unique nights.
 */ 
bool trackHasSufficientSupport(const std::vector<MopsDetection> &allDetections,
                               const Track &newTrack, const linkTrackletsConfig &searchConfig)
{
    if (newTrack.getComponentDetectionIndices().size() < 
        searchConfig.minDetectionsPerTrack) {
        return false;
    }
    /* count the number of NIGHTS of support.  Do this in a slightly
     * hackish way: if two successive detections are > .5 days apart
     * in time, then they must be from different nights.  So just
     * order all the detection MJDs and then walk across them looking
     * for unique nights.
     */
    const std::set<uint> trackDets = newTrack.getComponentDetectionIndices();
    std::set<double> allMjds;
    for (std::set<uint>::const_iterator detIter = trackDets.begin();
         detIter != trackDets.end(); detIter++) {
        allMjds.insert(allDetections.at(*detIter).getEpochMJD());
    }

    uint uniqueNightsSeen = 1;
    // this line could cause an exception if the program is run with
    // minDetectionsPerTrack <= 0 but that'd be totally insane.
    double lastDetTime = *(allMjds.begin());

    // note: we start at the second item, because we already looked at
    // the first. Ugly C++ syntax.
    std::set<double>::const_iterator mjdIter = allMjds.begin();
    for (++mjdIter; mjdIter != allMjds.end();  mjdIter++) {
        
        if (*mjdIter - lastDetTime > .5) {
            uniqueNightsSeen++;
        }
        lastDetTime = *mjdIter;            
    }
    
    if (uniqueNightsSeen < searchConfig.minUniqueNights) {
        return false;
    }
    return true;
}





/* 
 * The name no longer reflects exactly what this does.  It decides
 * whether the quality of the fit of the track model to the detections
 * is sufficiently good.   This is based on a combination of prob(chisq)
 * and the requirement of a physical range fit from the topocentric corrections.
 * Due to the definition of the fit function, a physical range value here is
 * NEGATIVE.  If a range was not fit, it will be zero.
 */
 
bool trackRmsIsSufficientlyLow(
    const std::vector<MopsDetection> &allDetections,
    const Track &newTrack, 
    const linkTrackletsConfig &searchConfig)
{

    bool ok = newTrack.getProbChisqRa() > searchConfig.trackMinProbChisq &&
        newTrack.getProbChisqDec() > searchConfig.trackMinProbChisq &&
        newTrack.getFitRange() <= 0.0;

    return ok;
}







/*
 * this is called when all endpoint nodes (i.e. model nodes) and
 * support nodes are leaves.  model nodes and support nodes are
 * expected to be mutually compatible.
 */
void buildTracksAddToResults(
    const std::vector<MopsDetection> &allDetections,
    const std::vector<Tracklet> &allTracklets,
    const linkTrackletsConfig &searchConfig,
    TreeNodeAndTime &firstEndpoint,
    TreeNodeAndTime &secondEndpoint,
    std::vector<TreeNodeAndTime> &supportNodes,
    TrackSet & results)
{

    if ((firstEndpoint.myTree->isLeaf() == false) ||
        (secondEndpoint.myTree->isLeaf() == false)) {
        LSST_EXCEPT(ProgrammerErrorException, 
            "buildTracksAddToResults got non-leaf nodes, must be a bug!");
    }
    for (uint i = 0; i < supportNodes.size(); i++) {
        if (supportNodes.at(i).myTree->isLeaf() == false) {
            LSST_EXCEPT(ProgrammerErrorException, 
           "buildTracksAddToResults got non-leaf nodes, must be a bug!"); 
        }
    }

    std::vector<PointAndValue<uint> >::const_iterator firstEndpointIter;
    std::vector<PointAndValue<uint> >::const_iterator secondEndpointIter;
    std::vector<TreeNodeAndTime>::const_iterator supportNodeIter;
    std::vector<PointAndValue <uint> >::const_iterator supportPointIter;
    for (firstEndpointIter = firstEndpoint.myTree->getMyData()->begin();
         firstEndpointIter != firstEndpoint.myTree->getMyData()->end();
         firstEndpointIter++) {

        for (secondEndpointIter = secondEndpoint.myTree->getMyData()->begin();
             secondEndpointIter != secondEndpoint.myTree->getMyData()->end();
             secondEndpointIter++) {
            
            /* figure out the rough quadratic track fitting the two
             * endpoints.  if error is too large, quit. Otherwise,
             * choose support points from the support nodes, using
             * best-fit first, and ignoring those too far off the
             * line.  If we get enough points, return a track.
            */

            // create a new track with these endpoints
            Track newTrack;
            
            uint firstEndpointTrackletIndex = firstEndpointIter->getValue();
            uint secondEndpointTrackletIndex = secondEndpointIter->getValue();
            

            
            newTrack.addTracklet(firstEndpointTrackletIndex, 
                                 allTracklets.at(firstEndpointTrackletIndex),
                                 allDetections);
            
            newTrack.addTracklet(secondEndpointTrackletIndex, 
                                 allTracklets.at(secondEndpointTrackletIndex),
                                 allDetections);
            // the 3 here says do NOT use the full form for ra and dec fit - use quadratic
            newTrack.calculateBestFitQuadratic(allDetections, 3);

            if (endpointTrackletsAreCompatible(allDetections, 
                                               newTrack,
                                               searchConfig)) {
                
                std::vector<uint> candidateTrackletIds;
                // put all support tracklet Ids in curSupportNodeData,
                // then call addDetectionsCloseToPredictedPositions
                for (supportNodeIter = supportNodes.begin(); 
                     supportNodeIter != supportNodes.end();
                     supportNodeIter++) {
                    const std::vector<PointAndValue <uint> > * 
                        curSupportNodeData;
                    if (!supportNodeIter->myTree->isLeaf()) {
                        throw LSST_EXCEPT(BadParameterException,
                                          std::string(__FUNCTION__) + 
                                          std::string(
                         ": received non-leaf node as support node."));
                    }
                    curSupportNodeData = supportNodeIter->myTree->getMyData(); 
                    for (supportPointIter  = curSupportNodeData->begin(); 
                         supportPointIter != curSupportNodeData->end();
                         supportPointIter++) {

                        candidateTrackletIds.push_back(
                            supportPointIter->getValue());
                    }
                }

                // Add support points if they are within thresholds.
                addDetectionsCloseToPredictedPositions(allDetections, 
                                                       allTracklets, 
                                                       candidateTrackletIds,
                                                       newTrack, 
                                                       searchConfig);

                
                // Final check to see if track meets requirements:
                // first sufficient support, then RMS fitting error

                if (trackHasSufficientSupport(allDetections, 
                                              newTrack,
                                              searchConfig)) {
                    // recalculate best-fit quadratic and check if RMS
                    // is sufficient
                    // the 'true' here says DO use the full form for ra fit if there are enough
                    // points to do so

                    newTrack.calculateBestFitQuadratic(allDetections, -1);
                    if (trackRmsIsSufficientlyLow(allDetections, 
                                                  newTrack, 
                                                  searchConfig)) {
#ifdef DEBUG
                        std::cout << "track passed rms" << std::endl;
#endif

                        results.insert(newTrack);

                    } else {
#ifdef DEBUG
                        std::cout << "track failed rms" << std::endl;
#endif
                    }
                } else {
#ifdef DEBUG
                    std::cout << "track has insufficient support" << std::endl;
#endif
                }

            }  else {
#ifdef DEBUG
                    std::cout << "endpoint tracklets incompatible" << std::endl;
#endif
                }  
    }
}


}















/* 
 * feb 17, 2011: update acc bounds using formulas reverse-engineered
 * from Kubica.
 */
bool updateAccBoundsReturnValidity(const TreeNodeAndTime &firstEndpoint, 
                                   const TreeNodeAndTime &secondEndpoint,
                                   double &aMinRa, double &aMaxRa, 
                                   double &aMinDec, double &aMaxDec)
{

    double dt = secondEndpoint.myTime.getMJD() - 
        firstEndpoint.myTime.getMJD();
    double dt2 = 2./(dt*dt);
    double dti = 1./(dt);

    TrackletTreeNode* A = firstEndpoint.myTree;
    TrackletTreeNode* B = secondEndpoint.myTree;

    double AmaxVRa, AminVRa, AmaxPRa, AminPRa;
    double BmaxVRa, BminVRa, BmaxPRa, BminPRa;
    double AmaxVDec, AminVDec, AmaxPDec, AminPDec;
    double BmaxVDec, BminVDec, BmaxPDec, BminPDec;

    // short-circuit right away if possible.
    if ((aMaxRa < aMinRa) || (aMaxDec < aMinDec)) {
        return false;
    }

    AmaxPRa = A->getUBounds()->at(POINT_RA);
    AminPRa = A->getLBounds()->at(POINT_RA);
    AmaxVRa = A->getUBounds()->at(POINT_RA_VELOCITY);
    AminVRa = A->getLBounds()->at(POINT_RA_VELOCITY);
    
    BmaxPRa = B->getUBounds()->at(POINT_RA);
    BminPRa = B->getLBounds()->at(POINT_RA);
    BmaxVRa = B->getUBounds()->at(POINT_RA_VELOCITY);
    BminVRa = B->getLBounds()->at(POINT_RA_VELOCITY);

    AmaxPDec = A->getUBounds()->at(POINT_DEC);
    AminPDec = A->getLBounds()->at(POINT_DEC);
    AmaxVDec = A->getUBounds()->at(POINT_DEC_VELOCITY);
    AminVDec = A->getLBounds()->at(POINT_DEC_VELOCITY);
    
    BmaxPDec = B->getUBounds()->at(POINT_DEC);
    BminPDec = B->getLBounds()->at(POINT_DEC);
    BmaxVDec = B->getUBounds()->at(POINT_DEC_VELOCITY);
    BminVDec = B->getLBounds()->at(POINT_DEC_VELOCITY);

    
    double tmpAcc;
    // max using vel test
    tmpAcc = (BmaxVRa - AminVRa) * dti;
    if (tmpAcc < aMaxRa) {
        aMaxRa = tmpAcc;
    }
    tmpAcc = (BmaxVDec - AminVDec) * dti;
    if (tmpAcc < aMaxDec) {
        aMaxDec = tmpAcc;
    }
    // min using velocity test
    tmpAcc = (BminVRa - AmaxVRa) * dti;
    if (tmpAcc > aMinRa) {
        aMinRa = tmpAcc;
    }
    tmpAcc = (BminVDec - AmaxVDec) * dti;
    if (tmpAcc > aMinDec) {
        aMinDec = tmpAcc;
    }
    // short circuit if possible
    if ((aMinRa > aMaxRa) || (aMinDec > aMaxDec)) {
        return false;
    }


    // max using pos/vel test 1
    tmpAcc = dt2 * (AmaxPRa - BminPRa + BmaxVRa * dt);
    if (tmpAcc < aMaxRa) {
        aMaxRa = tmpAcc;
    }
    tmpAcc = dt2 * (AmaxPDec - BminPDec + BmaxVDec * dt);
    if (tmpAcc < aMaxDec) {
        aMaxDec = tmpAcc;
    }
    // min using pos/vel test 1
    tmpAcc = dt2*(BminPRa - AmaxPRa - AmaxVRa * dt);
    if (tmpAcc > aMinRa) {
        aMinRa = tmpAcc;
    }
    tmpAcc = dt2*(BminPDec - AmaxPDec - AmaxVDec * dt);
    if (tmpAcc > aMinDec) {
        aMinDec = tmpAcc;
    }
    // short circuit if possible
    if ((aMinRa > aMaxRa) || (aMinDec > aMaxDec)) {
        return false;
    }


    // max using pos/vel test 2
    tmpAcc = dt2 * (BmaxPRa - AminPRa - AminVRa * dt);
    if (tmpAcc < aMaxRa) {
        aMaxRa = tmpAcc;
    }
    tmpAcc = dt2 * (BmaxPDec - AminPDec - AminVDec * dt);
    if (tmpAcc < aMaxDec) {
        aMaxDec = tmpAcc;
    }
    // min using pos/vel test 2
    tmpAcc = dt2*(AminPRa - BmaxPRa + BminVRa * dt);
    if (tmpAcc > aMinRa) {
        aMinRa = tmpAcc;
    }
    tmpAcc = dt2*(AminPDec - BmaxPDec + BminVDec * dt);
    if (tmpAcc > aMinDec) {
        aMinDec = tmpAcc;
    }
    // short circuit if possible
    if ((aMinRa > aMaxRa) || (aMinDec > aMaxDec)) {
        return false;
    }
    
    
    // we know maxAcc > minAcc because we didn't short-circuit above.
    return true;
}



/*
 * determine whether support node is compatible with endpoint nodes
 * using math borrowed from C linkTracklets.  Refines min,max
 * acceleration limits; we assume that these are then potentially used
 * for splitting child nodes of the support node in question
 */
bool areMutuallyCompatible(const TreeNodeAndTime &firstNode,
                           const TreeNodeAndTime &secondNode,
                           const TreeNodeAndTime &thirdNode,
                           const linkTrackletsConfig &searchConfig,
                           double &aMinRa, double &aMaxRa,
                           double &aMinDec, double &aMaxDec)
{

    bool firstPairCompat = 
        updateAccBoundsReturnValidity(firstNode, secondNode, 
                                      aMinRa, aMaxRa, aMinDec, aMaxDec);
    if (!firstPairCompat) 
        return false;

    bool secondPairCompat = 
        updateAccBoundsReturnValidity(secondNode, thirdNode, 
                                      aMinRa, aMaxRa, aMinDec, aMaxDec);
    return secondPairCompat;
}






bool areAllLeaves(const std::vector<TreeNodeAndTime> &nodeArray) {
    bool allLeaves = true;
    std::vector<TreeNodeAndTime>::const_iterator treeIter;
    uint count = 0;
    for (treeIter = nodeArray.begin(); 
         (treeIter != nodeArray.end() && (allLeaves == true));
         treeIter++) {
        if (treeIter->myTree->isLeaf() == false) {
            allLeaves = false;
        }
        count++;
    }
    return allLeaves;
}







bool supportTooWide(const TreeNodeAndTime& firstEndpoint, 
                    const TreeNodeAndTime& secondEndpoint,
                    const TreeNodeAndTime& supportNode) 
{
    /* odd-looking stuff with alpha based on test_and_add_support in
     * Kubica's linker.c.  the idea is to weight the expected size of a
     * node according to its proximity to either the first or last
     * endpoint . */
    double alpha = (supportNode.myTime.getMJD() - firstEndpoint.myTime.getMJD()) /
        (secondEndpoint.myTime.getMJD() - firstEndpoint.myTime.getMJD());
    
    // check the width of the support node in all 4 axes; compare with
    // width of
    for (uint i = 0; i < 4; i++) {
        double nodeWidth = supportNode.myTree->getUBounds()->at(i) - 
            supportNode.myTree->getLBounds()->at(i);
        double maxWidth = (1.0 - alpha) * 
            (firstEndpoint.myTree->getUBounds()->at(i) - 
             firstEndpoint.myTree->getLBounds()->at(i)) + 
            alpha * (secondEndpoint.myTree->getUBounds()->at(i) - 
                     secondEndpoint.myTree->getLBounds()->at(i));
        
        // this constant 4 is taken from Kubica's linker.c
        // test_and_add_support.  Beware the magic number!
        if (4.0 * nodeWidth > maxWidth) {
            return true;
        }
    }

    return false;
}







void splitSupportRecursively(const TreeNodeAndTime& firstEndpoint, 
                             const TreeNodeAndTime& secondEndpoint, 
                             bool requireLeaves,
                             const TreeNodeAndTime &supportNode, 
                             const linkTrackletsConfig &searchConfig, 
                             double accMinRa, double accMaxRa,
                             double accMinDec, double accMaxDec,
                             std::vector<TreeNodeAndTime> &newSupportNodes)
{

    if ((firstEndpoint.myTime.getMJD() >= supportNode.myTime.getMJD()) || 
        (supportNode.myTime.getMJD() > secondEndpoint.myTime.getMJD())) {
        throw LSST_EXCEPT(BadParameterException, "splitSupportRecursively got impossibly-ordered endpoints/support");
    }

    
    if (areMutuallyCompatible(firstEndpoint, supportNode,
                              secondEndpoint, searchConfig, 
                              accMinRa, accMaxRa,
                              accMinDec, accMaxDec)) {

        if (supportNode.myTree->isLeaf()) {
            newSupportNodes.push_back(supportNode);
        }

        else if (requireLeaves) {
            if (supportNode.myTree->hasLeftChild()) {
                TreeNodeAndTime leftTat(
                    supportNode.myTree->getLeftChild(), 
                    supportNode.myTime); 
                splitSupportRecursively(firstEndpoint, 
                                        secondEndpoint, 
                                        requireLeaves, 
                                        leftTat,
                                        searchConfig, 
                                        accMinRa, accMaxRa,
                                        accMinDec, accMaxDec,
                                        newSupportNodes);
            }
            if (supportNode.myTree->hasRightChild()) {
                TreeNodeAndTime rightTat(
                    supportNode.myTree->getRightChild(), 
                    supportNode.myTime);
                splitSupportRecursively(firstEndpoint, 
                                        secondEndpoint, 
                                        requireLeaves, 
                                        rightTat,
                                        searchConfig, 
                                        accMinRa, accMaxRa,
                                        accMinDec, accMaxDec,
                                        newSupportNodes);
            }
        }
        
        else {
            // we don't require leaves in output, but check to see if
            // we *should* split this node.  if not, add it to
            // output. Otherwise, recurse on its children.

            bool tooWide = supportTooWide(firstEndpoint, 
                                          secondEndpoint, 
                                          supportNode);
            
            if (tooWide) {
                if (supportNode.myTree->hasLeftChild()) {
                    TreeNodeAndTime leftTat(
                        supportNode.myTree->getLeftChild(), 
                        supportNode.myTime); 
                    splitSupportRecursively(firstEndpoint, 
                                            secondEndpoint, 
                                            requireLeaves, 
                                            leftTat,
                                            searchConfig, 
                                            accMinRa, accMaxRa,
                                            accMinDec, accMaxDec,
                                            newSupportNodes);
                }
                if (supportNode.myTree->hasRightChild()) {
                    TreeNodeAndTime rightTat(
                        supportNode.myTree->getRightChild(), 
                        supportNode.myTime);
                    splitSupportRecursively(firstEndpoint, 
                                            secondEndpoint, 
                                            requireLeaves, 
                                            rightTat,
                                            searchConfig, 
                                            accMinRa, accMaxRa,
                                            accMinDec, accMaxDec,
                                            newSupportNodes);
                }
            }
            else {
                newSupportNodes.push_back(supportNode);
            }
        }
    }
}





void filterAndSplitSupport(
    const TreeNodeAndTime& firstEndpoint, 
    const TreeNodeAndTime& secondEndpoint, 
    const std::vector<TreeNodeAndTime> &supportNodes, 
    const linkTrackletsConfig &searchConfig, 
    double accMinRa, double accMaxRa, 
    double accMinDec, double accMaxDec,
    std::vector<TreeNodeAndTime> &newSupportNodes) 
{

    // if the endpoints are leaves, require that we get all leaves in
    // the support nodes.
    bool endpointsAreLeaves = 
        firstEndpoint.myTree->isLeaf() && secondEndpoint.myTree->isLeaf();
    
    
    for (uint i = 0; i < supportNodes.size(); i++) {
        splitSupportRecursively(firstEndpoint, secondEndpoint, 
                                endpointsAreLeaves, 
                                supportNodes[i],
                                searchConfig, 
                                accMinRa, accMaxRa, accMinDec, accMaxDec,
                                newSupportNodes);
    }
    
    
}







double nodeWidth(TrackletTreeNode *node)
{
    double width = 1;
    for (uint i = 0; i < 4; i++) {
        width *= node->getUBounds()->at(i) - node->getLBounds()->at(i);
    }
    return width;    
}









unsigned int countImageTimes(const std::vector<TreeNodeAndTime> &nodes)
{
    std::set<unsigned int> imageTimes;
    std::vector<TreeNodeAndTime>::const_iterator nIter;
    for (nIter = nodes.begin(); nIter != nodes.end(); nIter++) {
        imageTimes.insert(nIter->myTime.getImageId());
    }
    return imageTimes.size();
}



/*
 * this is, roughly, the algorithm presented in http://arxiv.org/abs/astro-ph/0703475v1:
 * 
 * Efficient intra- and inter-night linking of asteroid detections
 * using kd-trees
 * 
 * The "proper" version of the algorithm (as presented in psuedocode
 * om the document) is a little different; basically, he has us
 * recurring only on endpoints, and splitting support nodes repeatedly
 * at each call, until they are no longer "too wide".  Unfortunately,
 * I found that neither my intuition nor Kubica's implementation made
 * clear the definition of "too wide".
 *
 * this implementation is more like the one Kubica did in his code. at
 * every step, we check all support nodes for compatibility, splitting
 * each one. we then split one model node and recurse.
 */
void doLinkingRecurse(const std::vector<MopsDetection> &allDetections,
                      const std::vector<Tracklet> &allTracklets,
                      const linkTrackletsConfig &searchConfig,
                      TreeNodeAndTime &firstEndpoint,
                      TreeNodeAndTime &secondEndpoint,
                      std::vector<TreeNodeAndTime> &supportNodes,
                      double accMinRa, double accMaxRa, 
                      double accMinDec, double accMaxDec,
                      TrackSet & results,
                      int iterationsTillSplit)
{

    firstEndpoint.myTree->addVisit();


    bool isValid = updateAccBoundsReturnValidity(firstEndpoint, 
                                                 secondEndpoint,
                                                 accMinRa, accMaxRa, 
                                                 accMinDec, accMaxDec); 

    if (isValid)
    {

        std::set<double> uniqueSupportMJDs;
        std::vector<TreeNodeAndTime> newSupportNodes;
        std::vector<TreeNodeAndTime>::iterator supportNodeIter;
        
        /* look through untested support nodes, find the ones that are
         * compatible with the model nodes, add their children to
         * newSupportNodes */

        if ((iterationsTillSplit <= 0) || 
            (firstEndpoint.myTree->isLeaf() && secondEndpoint.myTree->isLeaf())) {
            
            filterAndSplitSupport(firstEndpoint, secondEndpoint, 
                                  supportNodes, searchConfig, 
                                  accMinRa, accMaxRa, accMinDec, accMaxDec,
                                  newSupportNodes);
            iterationsTillSplit = ITERATIONS_PER_SPLIT;
        }
        else{
            newSupportNodes = supportNodes;
        }
        
        unsigned int nUniqueMJDs = countImageTimes(newSupportNodes);

        // we get at least 2 unique nights from endpoints, and 4
        // unique detections from endpoints.  add those in and see if
        // we have sufficient support.
        if (nUniqueMJDs + 2 >= searchConfig.minUniqueNights)
        {
            // we have enough model nodes, and enough support nodes.
            // if they are all leaves, then start building tracks.  if
            // they are not leaves, split one of them and recurse.

            if (firstEndpoint.myTree->isLeaf() && 
                secondEndpoint.myTree->isLeaf()) {
                
                buildTracksAddToResults(allDetections, 
                                        allTracklets, 
                                        searchConfig,
                                        firstEndpoint, 
                                        secondEndpoint, 
                                        newSupportNodes,
                                        results);
            }
            else {
                
                iterationsTillSplit -= 1;
                
                // find the "widest" node, where width is just the
                // product of RA range, Dec range, RA velocity range,
                // Dec velocity range.  we will split that node and
                // recurse.
                
                double firstEndpointWidth = nodeWidth(firstEndpoint.myTree);
                double secondEndpointWidth = nodeWidth(secondEndpoint.myTree);
                
                // don't consider splitting endpoint nodes which are
                // actually leaves! give them negative width to hack
                // selection process.
                if (firstEndpoint.myTree->isLeaf()) {
                    firstEndpointWidth = -1;
                }
                if (secondEndpoint.myTree->isLeaf()) {
                    secondEndpointWidth = -1;
                }
 
                // choose the widest model node, split it and recurse!
                if (firstEndpointWidth >= secondEndpointWidth) {

                    //"widest" node is first endpoint, recurse twice
                    // using its children in its place.
                    
                    if ((! firstEndpoint.myTree->hasLeftChild()) 
                        && (!firstEndpoint.myTree->hasRightChild())) {
                        throw LSST_EXCEPT(ProgrammerErrorException, 
     "Recursing in a leaf node (first endpoint), must be a bug!");
                    }

                    if (firstEndpoint.myTree->hasLeftChild())
                    {
                        TreeNodeAndTime newTAT(
                            firstEndpoint.myTree->getLeftChild(), 
                            firstEndpoint.myTime);
                        //doLinkingRecurseTime += timeSince(start);
                        //std::cout << "Recursing on left child of
                        //first endpoint." << std::endl;
                        doLinkingRecurse(allDetections, 
                                         allTracklets, 
                                         searchConfig,
                                         newTAT,secondEndpoint,
                                         newSupportNodes,
                                         accMinRa,
                                         accMaxRa,
                                         accMinDec,
                                         accMaxDec,
                                         results, 
                                         iterationsTillSplit); 
                        //std::cout << "Returned from recursion on
                        //left child of first endpoint." << std::endl;
                    }
                    
                    if (firstEndpoint.myTree->hasRightChild())
                    {
                        TreeNodeAndTime newTAT(
                            firstEndpoint.myTree->getRightChild(), 
                            firstEndpoint.myTime);
                        //std::cout << "recursing on right child of
                        //first endpoint.." << std::endl;
                        doLinkingRecurse(allDetections, 
                                         allTracklets, 
                                         searchConfig,
                                         newTAT,secondEndpoint,
                                         newSupportNodes,
                                         accMinRa,
                                         accMaxRa,
                                         accMinDec,
                                         accMaxDec,
                                         results, 
                                         iterationsTillSplit);  
                        //std::cout << "Returned from recursion on
                        //right child of first endpoint." << std::endl;
                    }
                }
                else {
                    //"widest" node is second endpoint, recurse twice
                    // using its children in its place
                    
                    if ((!secondEndpoint.myTree->hasLeftChild()) 
                        && (!secondEndpoint.myTree->hasRightChild())) {
                        throw LSST_EXCEPT(ProgrammerErrorException, 
            "Recursing in a leaf node (second endpoint), must be a bug!");
                    }

                    if (secondEndpoint.myTree->hasLeftChild())
                    {
                        TreeNodeAndTime newTAT(
                            secondEndpoint.myTree->getLeftChild(), 
                            secondEndpoint.myTime);
                        //std::cout << "Recursing on left child of
                        //second endpoint." << std::endl;
                        doLinkingRecurse(allDetections, 
                                         allTracklets, 
                                         searchConfig,
                                         firstEndpoint,
                                         newTAT,
                                         newSupportNodes,
                                         accMinRa,
                                         accMaxRa,
                                         accMinDec,
                                         accMaxDec,
                                         results, 
                                         iterationsTillSplit);
                        //std::cout << "Returned from recursion on
                        //left child of second endpoint." << std::endl;
                    }
                    
                    if (secondEndpoint.myTree->hasRightChild())
                    {
                        TreeNodeAndTime newTAT(
                            secondEndpoint.myTree->getRightChild(), 
                            secondEndpoint.myTime);
                        //std::cout << "Recursing on right child of
                        //second endpoint." << std::endl;
                        doLinkingRecurse(allDetections, 
                                         allTracklets, 
                                         searchConfig,
                                         firstEndpoint,
                                         newTAT,
                                         newSupportNodes,
                                         accMinRa,
                                         accMaxRa,
                                         accMinDec,
                                         accMaxDec,
                                         results, 
                                         iterationsTillSplit);
                        //std::cout << "Returned from recursion on
                        //right child of second endpoint." << std::endl;
                        
                    }
                }
            }                        
        }
    }
}



class WorkItem {
public:
    std::map<ImageTime, TrackletTree>::const_iterator firstEndpoint;
    std::map<ImageTime, TrackletTree>::const_iterator secondEndpoint;    
};


void doLinking(const std::vector<MopsDetection> &allDetections,
               std::vector<Tracklet> &allTracklets,
               const linkTrackletsConfig &searchConfig,
               std::map<ImageTime, TrackletTree > &trackletTimeToTreeMap,
               TrackSet &results)
{
    /* for every pair of trees, using the set of every intermediate
     * (temporally) tree as a set possible support nodes, call the
     * recursive linker.
     */ 
    /*
     * implementation detail: note that there is a difference between
     * a KDTree` and a KDTreeNode.  KDTrees are the "outer layer" and
     * don't really hold data; KDTreeNodes are where the data lives.
     * 
     * In other programs, KDTreeNodes are hidden from the user, but
     * linkTracklets will use them directly.
     *
     * In doLinking_recurse, we will only deal with KDTreeNodes. So in
     * doLinking, we start the searching with the head (i.e. root)
     * node of each tree.
     */
    bool DEBUG = false;
    unsigned int imagePairs = 0;

    unsigned int numImages = trackletTimeToTreeMap.size();

    std::map<ImageTime, TrackletTree >::const_iterator firstEndpointIter;

    /* OMP doesn't deal well with for loops on iterators. Need to
       create an array of work items to do and loops on that. */
    std::vector<WorkItem> allWork;
    

    for (firstEndpointIter = trackletTimeToTreeMap.begin(); 
         firstEndpointIter != trackletTimeToTreeMap.end(); 
         firstEndpointIter++)
    {
        std::map<ImageTime, TrackletTree >::const_iterator 
            secondEndpointIter;
        std::map<ImageTime, TrackletTree >::const_iterator 
            afterFirstIter = firstEndpointIter;
        afterFirstIter++;


        /* check if the user wants us to look for tracks starting at this time */
        if (( !searchConfig.restrictTrackStartTimes ) || 
            (firstEndpointIter->first.getMJD() 
             <= searchConfig.latestFirstEndpointTime)) {

            for (secondEndpointIter = afterFirstIter; 
                 secondEndpointIter != trackletTimeToTreeMap.end(); 
                 secondEndpointIter++)
            {
                
                /* check if the user wants us to look for tracks
                 * ending at this time */
                if ((!searchConfig.restrictTrackEndTimes) || 
                    (secondEndpointIter->first.getMJD() 
                     >= searchConfig.earliestLastEndpointTime)) {
                
                    /* if there is sufficient time between the first
                       and second nodes, then try doing the recursive
                       linking using intermediate times as support
                       nodes.
                    */
                    
                    if (secondEndpointIter->first.getMJD() 
                        - firstEndpointIter->first.getMJD() 
                        >= searchConfig.minEndpointTimeSeparation) {

                        // get all intermediate points as support
                        // nodes.
                
                        /* note that std::maps are sorted by their
                         * key, which in this case is time.  ergo
                         * between firstEndpointIter and
                         * secondEndpointIter is *EVERY* tree (and
                         * ergo every tracklet) which happened between
                         * the first endpoint's tracklets and the
                         * second endpoint's tracklets.
                         */
                
                        uint nSupport = 0;
                        std::map<ImageTime, TrackletTree >::const_iterator 
                            supportPointIter;
                        // look for support; break if we find any at all.
                        for (supportPointIter = afterFirstIter;
                             supportPointIter != secondEndpointIter;
                             supportPointIter++) {
                    
                            /* don't pass along second tracklets which
                               are 'too close' to the endpoints; see
                               linkTracklets.h for more comments */
                            double firstToSup = supportPointIter->first.getMJD() 
                                - firstEndpointIter->first.getMJD(); 
                            double supToSecond =  secondEndpointIter->first.getMJD() 
                                - supportPointIter->first.getMJD();
                            if ((firstToSup > 
                                 searchConfig.minSupportToEndpointTimeSeparation) 
                                && 
                                (supToSecond > 
                                 searchConfig.minSupportToEndpointTimeSeparation)) 
                            {
                                nSupport += 1;
                            }
                        }
                
                        if (nSupport > 0) {
                            // check to see if the top-level tree
                            // nodes are really compatible.  we can
                            // probably cut out a lot of overhead if
                            // we only allocate hard work, and no
                            // trivial work.

                            imagePairs++;
                            TreeNodeAndTime firstEndpoint(firstEndpointIter->second.getRootNode(), 
                                                          firstEndpointIter->first);
                            TreeNodeAndTime secondEndpoint(secondEndpointIter->second.getRootNode(),
                                                           secondEndpointIter->first);
                            double accMaxRa = searchConfig.maxRAAccel;
                            double accMinRa = accMaxRa * -1;
                            double accMaxDec = searchConfig.maxDecAccel;
                            double accMinDec = accMaxDec * -1;
                            bool isValid = updateAccBoundsReturnValidity(firstEndpoint, 
                                                                         secondEndpoint,
                                                                         accMinRa, 
                                                                         accMaxRa, 
                                                                         accMinDec, 
                                                                         accMaxDec); 

                            if (isValid)
                                {
                                    // create a work item here and
                                    // pass it off.
                                    WorkItem newWork;
                                    newWork.firstEndpoint = firstEndpointIter;
                                    newWork.secondEndpoint = secondEndpointIter;
                                    allWork.push_back(newWork); 
                                }
                            else {
                                TrackSet tmpRes;
                                // we should get NO results. if we do then panic.
        std::vector<TreeNodeAndTime > supportPoints;
        std::map<ImageTime, TrackletTree >::const_iterator 
            supportPointIter;
        std::map<ImageTime, TrackletTree>::const_iterator afterFirstIter;
        afterFirstIter = firstEndpointIter;
        afterFirstIter++;
        for (supportPointIter = afterFirstIter;
             supportPointIter != secondEndpointIter;
             supportPointIter++) {
            
            /* don't pass along second tracklets which
               are 'too close' to the endpoints; see
               linkTracklets.h for more comments */
            double firstToSup = supportPointIter->first.getMJD() 
                - firstEndpointIter->first.getMJD(); 
            double supToSecond =  secondEndpointIter->first.getMJD() 
                - supportPointIter->first.getMJD();
            if ((firstToSup > 
                 searchConfig.minSupportToEndpointTimeSeparation) 
                && 
                (supToSecond > 
                 searchConfig.minSupportToEndpointTimeSeparation)) 
            {
                
                TreeNodeAndTime tmpTAT(supportPointIter->second.getRootNode(),
                                       supportPointIter->first);
                supportPoints.push_back(tmpTAT);
            }
        }
        // use tid to choose an output vector for just this thread to use.
        int tid;
        tid = omp_get_thread_num();

        doLinkingRecurse(allDetections,
                         allTracklets, 
                         searchConfig,
                         firstEndpoint, 
                         secondEndpoint,
                         supportPoints,  
                         searchConfig.maxRAAccel*-1.,
                         searchConfig.maxRAAccel,
                         searchConfig.maxDecAccel*-1.,
                         searchConfig.maxDecAccel,
                         tmpRes, 
                         ITERATIONS_PER_SPLIT);
        
        if (tmpRes.size() != 0) {
            std::cout << "WTF?! endpoints are not compatible but found " << tmpRes.size() << " tracks?!" << std::endl;

            accMaxRa = searchConfig.maxRAAccel;
            accMinRa = accMaxRa * -1;
            accMaxDec = searchConfig.maxDecAccel;
            accMinDec = accMaxDec * -1;

            bool isValid = updateAccBoundsReturnValidity(firstEndpoint, 
                                                         secondEndpoint,
                                                         accMinRa, accMaxRa, 
                                                         accMinDec, accMaxDec); 
            // call it again so the debugger can watch.
            doLinkingRecurse(allDetections,
                             allTracklets, 
                             searchConfig,
                             firstEndpoint, 
                             secondEndpoint,
                             supportPoints,  
                             searchConfig.maxRAAccel*-1.,
                             searchConfig.maxRAAccel,
                             searchConfig.maxDecAccel*-1.,
                             searchConfig.maxDecAccel,
                             tmpRes, 
                             ITERATIONS_PER_SPLIT);
        
        }
                            }
                        }
                    }
                }
            }
        }
    }

    std::vector<TrackSet*> localResults;
    
#pragma omp parallel shared(localResults)
    {
        int nthreads, tid;
        nthreads = omp_get_num_threads();
        tid = omp_get_thread_num();

        if(tid == 0) {

            std::cout << "Number of threads " << nthreads << std::endl;
            std::cout << "Number of endpoint pairs: " << imagePairs << std::endl;
            std::cout << "Number of post-filtered work items: " << allWork.size() << std::endl;

            // main thread sets up localResults vector for other threads to use
            for (int i = 0; i < nthreads; i++) {
                // configure each local result set to hold tracks in
                // memory (forever, not a temporary cache) and we will
                // explicitly flush to disk at the end.
                TrackSet *tmpTS = new TrackSet(searchConfig.outputFile, false, 0);
                localResults.push_back(tmpTS);
            }
        }
    }

#pragma omp parallel for schedule(dynamic, 128)
    for (uint i = 0; i < allWork.size(); i++) {
        std::map<ImageTime, TrackletTree>::const_iterator 
            firstEndpointIter, secondEndpointIter;
        firstEndpointIter = allWork[i].firstEndpoint;
        secondEndpointIter = allWork[i].secondEndpoint;
        
        TreeNodeAndTime firstEndpoint(firstEndpointIter->second.getRootNode(), 
                                      firstEndpointIter->first);
        TreeNodeAndTime secondEndpoint(secondEndpointIter->second.getRootNode(),
                                       secondEndpointIter->first);        

        std::vector<TreeNodeAndTime > supportPoints;
        std::map<ImageTime, TrackletTree >::const_iterator 
            supportPointIter;
        std::map<ImageTime, TrackletTree>::const_iterator afterFirstIter;
        afterFirstIter = firstEndpointIter;
        afterFirstIter++;
        for (supportPointIter = afterFirstIter;
             supportPointIter != secondEndpointIter;
             supportPointIter++) {
            
            /* don't pass along second tracklets which
               are 'too close' to the endpoints; see
               linkTracklets.h for more comments */
            double firstToSup = supportPointIter->first.getMJD() 
                - firstEndpointIter->first.getMJD(); 
            double supToSecond =  secondEndpointIter->first.getMJD() 
                - supportPointIter->first.getMJD();
            if ((firstToSup > 
                 searchConfig.minSupportToEndpointTimeSeparation) 
                && 
                (supToSecond > 
                 searchConfig.minSupportToEndpointTimeSeparation)) 
            {
                
                TreeNodeAndTime tmpTAT(supportPointIter->second.getRootNode(),
                                       supportPointIter->first);
                supportPoints.push_back(tmpTAT);
            }
        }
        // use tid to choose an output vector for just this thread to use.
        int tid;
        tid = omp_get_thread_num();

        doLinkingRecurse(allDetections,
                         allTracklets, 
                         searchConfig,
                         firstEndpoint, 
                         secondEndpoint,
                         supportPoints,  
                         searchConfig.maxRAAccel*-1.,
                         searchConfig.maxRAAccel,
                         searchConfig.maxDecAccel*-1.,
                         searchConfig.maxDecAccel,
                         *(localResults[tid]), 
                         ITERATIONS_PER_SPLIT);
        
    }

#pragma omp parallel shared(localResults)
    {
        int tid;
        tid = omp_get_thread_num();

        if(tid == 0) {
            
            std::cout << "Master thread: collecting results " << std::endl;

            // main thread sets up localResults vector for other threads to use
            for (int i = 0; i < localResults.size(); i++) {
                // the destructor will purge to file and clear the memory.
                delete localResults[i];
            }
        }
    }
    
}








TrackSet* linkTracklets(std::vector<MopsDetection> &allDetections,
                        std::vector<Tracklet> &queryTracklets,
                        const linkTrackletsConfig &searchConfig) {
    TrackSet * toRet;
    if (searchConfig.outputMethod == RETURN_TRACKS) {
        toRet = new TrackSet();
    }
    else if (searchConfig.outputMethod == IDS_FILE) {
        toRet = new TrackSet(searchConfig.outputFile, 
                             false, 
                             0);
    }
    else if (searchConfig.outputMethod == IDS_FILE_WITH_CACHE) {
        toRet = new TrackSet(searchConfig.outputFile, 
                             true, 
                             searchConfig.outputBufferSize);
    }
    else {
        throw LSST_EXCEPT(BadParameterException, 
      "linkTracklets: got unknown or unimplemented output method.");
    }


    /*create a sorted list of KDtrees, each tree holding tracklets
      with unique start times (times of first detection in the
      tracklet).
      
      the points in the trees are in [RA, Dec, RAVelocity,
      DecVelocity] and the returned keys are indices into
      queryTracklets.
    */
    if (searchConfig.myVerbosity.printStatus) {
        std::cout << "Recentering all detections on (180, 0)." << std::endl;
    }
    recenterDetections(allDetections, searchConfig);
    
    if (searchConfig.myVerbosity.printStatus) {
        std::cout << "Setting tracklet velocities." << std::endl;
    }
    setTrackletVelocities(allDetections, queryTracklets);

    if (searchConfig.myVerbosity.printStatus) {
        std::cout << "Sorting tracklets by image time and creating trees." << std::endl;
    }
    std::map<ImageTime, TrackletTree > trackletTimeToTreeMap;    
    makeTrackletTimeToTreeMap(allDetections, 
                              queryTracklets, 
                              trackletTimeToTreeMap, 
                              searchConfig);
    if (searchConfig.myVerbosity.printStatus) {
        std::cout << "Doing the linking." << std::endl;
    }
    doLinking(allDetections, 
              queryTracklets, 
              searchConfig, 
              trackletTimeToTreeMap, 
              *toRet);
    if (searchConfig.myVerbosity.printStatus) {
        std::cout << "Finished linking." << std::endl;
    }
    
    
    return toRet;
}


void calculateTopoCorr(std::vector<MopsDetection> &allDetections,
                       const linkTrackletsConfig &searchConfig) {

    std::vector<MopsDetection>::iterator detIter;
    
    for (detIter = allDetections.begin();
         detIter != allDetections.end();
         detIter++) {
        detIter->calculateTopoCorr();
    }


}



}} //close lsst::mops

// -*- LSST-C++ -*-
/* jonathan myers */

// time headers needed for benchmarking performance
#include <ctime>
#include <iomanip>
#include <map>
#include <gsl/gsl_multifit.h>
// even more GSL needed for new rotation matrix stuff
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <time.h>
#include <algorithm>

#include "lsst/mops/rmsLineFit.h"
#include "lsst/mops/daymops/linkTracklets/linkTracklets.h"
#include "lsst/mops/Exceptions.h"
#include "lsst/mops/KDTree.h"
#include "lsst/mops/daymops/linkTracklets/TrackletTree.h"




/* taking a queue from Kubica, it's only once per ITERATIONS_PER_SPLIT
 * calls to doLinkingRecurse that we actually split the (non-leaf)
 * support nodes. The idea is to avoid redundantly calculating whether
 * a set of support nodes is compatible.
 */
#define ITERATIONS_PER_SPLIT 0

/* 
the following flags, if set to 'true, will enable some debugging
checks which use brute-force searching and ground truth data to alert
the user if an object is going to be missed.  (not foolproof, but
pretty good).  Will make searching INSANELY slow.
 */
#define CHEAT_AND_DO_CORRECTNESS_CHECKS_AT_SUPPORT_POINTS false
#define CHEAT_AND_DO_CORRECTNESS_CHECKS_BY_ENDPOINT false
#define CHEAT_AND_REPORT_FINDABLE_OBJECTS false
#define CHEAT_AND_DO_TOP_LEVEL_CHECK false

#define RACE_TO_MAX_COMPATIBLE false
#define MAX_COMPATIBLE_TO_FIND 500000


#define POINT_RA           0
#define POINT_DEC          1
#define POINT_RA_VELOCITY   2
#define POINT_DEC_VELOCITY  3

#define uint unsigned int 

namespace lsst {
    namespace mops {



// Globals for measuring runtime.  Kill these off eventually!

double modifyWithAccelerationTime, 
    positionAndVelocityRangesOverlapAfterAccelerationTime, 
    calculateBestFitQuadraticTime,
    trackMeetsRequirementsTime, 
    buildTracksAddToResultsTime, 
    areAllLeavesTime, 
    nodeWidthTime, 
    addAllDetectedObjectsToSetTime, 
    doLinkingRecurseTime;

int doLinkingRecurseVisits, 
    buildTracksAddToResultsVisits, 
    compatibleEndpointsFound;
int rejectedOnVelocity, 
    rejectedOnPosition, 
    wereCompatible, 
    rejectedOnLackOfSupport;

int cacheHits, 
    cacheMisses;








/* **********************************************************************

* A FEW SIMPLE, LOCAL CLASSES: mostly for readability and to allow the
  cache of KDTree Node, Image time, second image Time -> bounding box

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


/*
 * this class is only used as a key to the cache of TreeNode, img time
 * T -> node bounding box at time T.
 *
 * note that KDTreeNodes have unique IDs *WITHIN THEIR TREE* but every
 * tree will contain a node with ID 1.  Tree node ID + Image Time
 * yields a unique identity.
 *
 */
class GlobalNodeIdAndProjectionId {
public:
    GlobalNodeIdAndProjectionId(uint _treeNodeId, 
                                uint _treeImgTime, 
                                uint _projectionTime) {
        treeNodeId = _treeNodeId;
        imgId = _treeImgTime;
        projectionId = _projectionTime;
    }
    GlobalNodeIdAndProjectionId() {
        treeNodeId = 0;
        imgId = 0;
    }

    void setImageId(uint newId) {
        imgId = newId;
    }

    void setTreeNodeId(uint newId) {
        treeNodeId = newId;
    }
    void setProjectionId(uint newId) {
        projectionId = newId;
    }

    uint getImageId() const {
        return imgId;
    }
    uint getTreeNodeId() const {
        return treeNodeId;
    }
    uint getProjectionId() const {
        return projectionId;
    }
    
    bool operator== (const GlobalNodeIdAndProjectionId &other) const {
        return (other.getImageId() == imgId) && 
            (other.getTreeNodeId() == treeNodeId) &&
            (other.getProjectionId() == projectionId);
    }

    bool operator< (const GlobalNodeIdAndProjectionId &other) const {
        // sort first on start image ID, then on
        // tree node ID, then on projection time ID.

        if (imgId == other.getImageId()) {

            if (treeNodeId == other.getTreeNodeId()) {

                return (projectionId < other.getProjectionId());
            }
            else {
                return (treeNodeId < other.getTreeNodeId());
            }
        }
        else {
            return (imgId < other.getImageId());
        }
    }

private:
    uint imgId;
    uint treeNodeId;
    uint projectionId;
};








class TreeNodeAndTime {
public:
    TreeNodeAndTime(TrackletTreeNode * tree, ImageTime i) {
        myTree = tree;
        myTime = i;
    }
    TrackletTreeNode * myTree;
    ImageTime myTime;

};










/* *****************************************************
* These functions are for debugging/ diagnostics only.  
********************************************************/


// for debugging and calculating timing info
double getTimeElapsed(clock_t priorEvent)
{
     return ( std::clock() - priorEvent ) / (double)CLOCKS_PER_SEC;
}
double timeSince(clock_t priorEvent)
{
     return ( std::clock() - priorEvent ) / (double)CLOCKS_PER_SEC;
}

void debugPrintTimingInfo(const TrackSet &results);



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










void debugPrintTimingInfo(const TrackSet &results)
{

    std::cout << " so far we've found " << results.size() << " tracks.\n";
    std::cout << "TIMING STATS: \n---------------------\n";
    std::cout << "modifyWithAcceleration:\t" << modifyWithAccelerationTime << "sec\n"; 
    std::cout << "positionAndVelocityRangesOverlapAfterAcceleration:\t" << positionAndVelocityRangesOverlapAfterAccelerationTime << "sec\n"; 
    std::cout << "trackMeetsRequirements:\t" << trackMeetsRequirementsTime << "sec\n"; 
    std::cout << "buildTracksAddToResults:\t" << buildTracksAddToResultsTime << "sec\n"; 
    std::cout << "areAllLeaves:\t" << areAllLeavesTime << "sec\n"; 
    std::cout << "nodeWidth:\t" << nodeWidthTime << "sec\n"; 
    std::cout << "addAllDetectedObjectsToSet:\t" << addAllDetectedObjectsToSetTime << "sec\n"; 
    std::cout << "doLinkingRecurse:\t" << doLinkingRecurseTime << "sec\n";
                                
    std::cout << "\n\nVisited doLinkingRecurse " << doLinkingRecurseVisits  << " times.\n";
    std::cout << "Visited buildTracksAddToResults " << buildTracksAddToResultsVisits << " times.\n";
    std::cout << "found " << compatibleEndpointsFound << " compatible endpoint pairs.\n";
    std::cout << "generated " << results.size() << " Tracks.\n";
                                
    std::cout << "cache had " << cacheHits << " hits.\n";
    std::cout << "cache had " << cacheMisses << " misses.\n";
                                
    std::cout << "While examining tracklets, found " << wereCompatible 
              << " compatibilities, counted separately in RA and Dec.\n";
    std::cout << "   - rejected " << rejectedOnVelocity << " based on velocity.\n";
    std::cout << "   - rejected " << rejectedOnPosition 
              << " based on position, after finding the compatible in velocity.\n";

    std::cout << "In recursion, terminated " << rejectedOnLackOfSupport 
              << " times due to insufficient support points.\n";

}



void initDebugTimingInfo()
{

    rejectedOnLackOfSupport = 0;
    rejectedOnVelocity = 0;
    rejectedOnPosition = 0;
    wereCompatible = 0;
    cacheHits = 0; 
    cacheMisses = 0;
    doLinkingRecurseVisits = 0;
    buildTracksAddToResultsVisits = 0;
    compatibleEndpointsFound = 0;
    modifyWithAccelerationTime = 0; 
    positionAndVelocityRangesOverlapAfterAccelerationTime = 0; 
    calculateBestFitQuadraticTime = 0;
    trackMeetsRequirementsTime = 0; 
    buildTracksAddToResultsTime = 0; 
    areAllLeavesTime = 0; 
    nodeWidthTime = 0; 
    addAllDetectedObjectsToSetTime = 0; 
    doLinkingRecurseTime = 0;    
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
    /* first, make sure that all data falls along a contiguous
     * region.  this will probably cause downstream failures in
     * subtle, hard-to-find ways if we start using full-sky data,
     * because our accBounds and areMutuallyCompatible math
     * doesn't deal with wrap-around. Please don't use this code
     * after DC3. */
    
    double firstRa, firstDec;
    firstRa = allDetections[0].getRA();
    firstDec = allDetections[0].getDec();
    double minRa = firstRa;
    double maxRa = firstRa;
    double minDec = firstDec;
    double maxDec = firstDec;
    for (unsigned int i = 1; i < allDetections.size(); i++) {
        double thisRa, thisDec;
        thisRa  = allDetections[i].getRA();
        thisDec = allDetections[i].getDec();
        while (firstRa - thisRa > 180.) {
            thisRa += 360.;
        }
        while (firstRa - thisRa < -180.) {
            thisRa -= 360.;
        }
        while (firstDec - thisDec > 180.) {
            thisDec += 360.;
        }
        while (firstDec - thisDec < -180.) {
            thisDec -= 360.;
        }
        allDetections[i].setRA(thisRa);
        allDetections[i].setDec(thisDec);
        
        minRa = minOfTwo(thisRa, minRa);
        maxRa = maxOfTwo(thisRa, maxRa);

        minDec = minOfTwo(thisDec, minDec);
        maxDec = maxOfTwo(thisDec, maxDec);
    }

    if ((maxRa - minRa >= 180.) || 
        (maxDec - minDec >= 180.)) {
        LSST_EXCEPT(KnownShortcomingException,
                    "Detections do not fall on contiguous (180,180) degree range. Math is known to fail in this case.");
    }

    if (searchConfig.myVerbosity.printBoundsInfo) {
        std::cout << "   Bounds were (in deg) R=[" << 
            minRa << ",  " << maxRa << "], D=[ " <<
            minDec << ",  " << maxDec << "]" << std::endl;
    }
    
    /* now use a rotation matrix to rotate all data to center on (180, 0). */
    Constants c;
    double raOld = c.deg_to_rad()*((minRa + maxRa) / 2.);
    double decOld = c.deg_to_rad()*((minDec + maxDec) / 2.);
    double newCenterRa = c.deg_to_rad()*180.;
    double newCenterDec = c.deg_to_rad()*0.;
    double dDec = newCenterDec - decOld;
    double tmp1[] = { 0., 0., 0., 
                      0., 0., 0.,
                      0., 0., 1.};
    tmp1[0] =  cos(raOld);
    tmp1[1] =  sin(raOld);
    tmp1[3] = -sin(raOld);
    tmp1[4] =  cos(raOld);
    double tmp2[] = {0., 0., 0.,
                     0., 1., 0.,
                     0., 0., 0.};
    tmp2[0] =  cos(dDec);
    tmp2[2] = -sin(dDec);
    tmp2[6] =  sin(dDec);
    tmp2[8] =  cos(dDec);
    double tmp3[] = {0., 0., 0.,
                     0., 0., 0.,
                     0., 0., 1.};
    tmp3[0] =  cos(-newCenterRa);
    tmp3[1] =  sin(-newCenterRa);
    tmp3[3] = -sin(-newCenterRa);
    tmp3[4] =  cos(-newCenterRa);

    double partial[] = {0., 0., 0.,
                        0., 0., 0.,
                        0., 0., 0.};

    gsl_matrix_view gTmp1 = gsl_matrix_view_array(tmp1, 3, 3);
    gsl_matrix_view gTmp2 = gsl_matrix_view_array(tmp2, 3, 3);
    gsl_matrix_view gTmp3 = gsl_matrix_view_array(tmp3, 3, 3);
    gsl_matrix_view gPart = gsl_matrix_view_array(partial, 3, 3);
    gsl_matrix * gRot  = gsl_matrix_alloc(3, 3);

    int err;

    // gPart = gTmp3 . gTmp2
    err = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                          1.0, &gTmp3.matrix, &gTmp2.matrix,
                          0.0, &gPart.matrix);
    if (err != 0) {
        LSST_EXCEPT(GSLException, "Error from gsl_blas_dgemm!");
    }
    // gRot = gPar . gTmp1
    err = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                          1.0, &gPart.matrix, &gTmp1.matrix,
                          0.0, gRot);
    if (err != 0) {
        LSST_EXCEPT(GSLException, "Error from gsl_blas_dgemm!");
    }

    // gRot is now our rotation matrix.

    gsl_vector *p0 = gsl_vector_alloc(3);
    gsl_vector *pR = gsl_vector_alloc(3);
    if ((p0 == NULL) || (pR == NULL)) {
        LSST_EXCEPT(GSLException, "Got NULL from gsl_vector_alloc!");
    }
    
    // reset center Ra, center Dec to degrees; 
    // we will use them to choose the numbering of RAs and Decs.
    newCenterRa = c.rad_to_deg()*newCenterRa;
    newCenterDec = c.rad_to_deg()*newCenterDec;

    // do the actual rotation; convert to cartesian, pR = gRot . p0,
    // convert from cartesian to RA/Dec again.
    for (uint i = 0; i < allDetections.size(); i++){
        double x,y,z;
        toCartesian_deg(allDetections[i].getRA(),
                        allDetections[i].getDec(),
                        x, y, z);
        gsl_vector_set(p0, 0, x);
        gsl_vector_set(p0, 1, y);
        gsl_vector_set(p0, 2, z);
        err = gsl_blas_dgemv(CblasNoTrans, 1.0, gRot, p0, 0., pR);
        if (err != 0) {
            LSST_EXCEPT(GSLException, "Error from gsl_blas_dgemm!\n");
        }
        x = gsl_vector_get(pR, 0);
        y = gsl_vector_get(pR, 1);
        z = gsl_vector_get(pR, 2);
        double newRa, newDec;
        toRaDec_deg(x, y, z, newRa, newDec);

        while (newCenterRa - newRa > 180.) {
            newRa += 360.;
        }
        while (newCenterRa - newRa < -180.) {
            newRa -= 360.;
        }
        while (newCenterDec - newDec > 180.) {
            newDec += 360.;
        }
        while (newCenterDec - newDec < -180.) {
            newDec -= 360.;
        }

        allDetections[i].setRA(newRa);
        allDetections[i].setDec(newDec);

        if (i == 0) {
            minRa = newRa;
            maxRa = newRa;
            minDec = newDec;
            maxDec = newDec;
        }

        minRa = minOfTwo(minRa, newRa);
        maxRa = maxOfTwo(maxRa, newRa);
        minDec = minOfTwo(minDec, newDec);
        maxDec = maxOfTwo(maxDec, newDec);
    }
    gsl_vector_free(p0);
    gsl_vector_free(pR);
    gsl_matrix_free(gRot);
    
    if (searchConfig.myVerbosity.printBoundsInfo) {
        std::cout << "   Bounds are (in deg) R=[" << 
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
                                        DecSlopeAndOffset);
        
        curTracklet->velocityRA = RASlopeAndOffset.at(0);
        curTracklet->velocityDec = DecSlopeAndOffset.at(0);
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
                  << " tracklets.\n";

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
                  << " image times.\n";
    
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
 * update position and velocity given acceleration over time.
 */
void modifyWithAcceleration(double &position, double &velocity, 
                            double acceleration, double time)
{
    double start = std::clock();
    // use good ol' displacement = vt + .5a(t^2) 
    double newPosition = position + velocity*time 
        + .5*acceleration*(time*time);
    double newVelocity = velocity + acceleration*time;
    position = newPosition;
    velocity = newVelocity;
    modifyWithAccelerationTime += timeSince(start);
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
                           double &minR, double &maxR,
                           double &minD, double &maxD)
{
    bool groundTruthShouldKeep = false;
    if ((firstNode.myTree->hasTracklet(79)) && (thirdNode.myTree->hasTracklet(116730))) {
        if (secondNode.myTree->hasTracklet(57930)) {
            groundTruthShouldKeep = true;
        }
    }

    for (uint whichPair = 0; whichPair < 2; whichPair++) {
        const TrackletTreeNode * A;
        const TrackletTreeNode * B;
        double aTime;
        double bTime;

        if (whichPair == 0) {
            A = firstNode.myTree;
            aTime = firstNode.myTime.getMJD();
            B = secondNode.myTree;
            bTime = secondNode.myTime.getMJD();            
        }
        else {
            A = secondNode.myTree;
            aTime = secondNode.myTime.getMJD();
            B = thirdNode.myTree;
            bTime = thirdNode.myTime.getMJD();
        }
    

        double dt = bTime - aTime;
        double dt2 = 2.0 / (dt * dt);
        double dti = 1.0/dt;
        
        for (uint axis = 0; axis < 2; axis++) {
            double parentMin, parentMax;
            double AmaxV, AminV, AmaxP, AminP;
            double BmaxV, BminV, BmaxP, BminP;
            unsigned int pos, vel;
            if (axis == 0) {
                pos=POINT_RA;
                vel=POINT_RA_VELOCITY;
                parentMin = minR;
                parentMax = maxR;
            }
            else {
                pos=POINT_DEC;
                vel=POINT_DEC_VELOCITY;
                parentMin = minD;
                parentMax = maxD;
            }

            // short circuit ASAP if this won't work
            if (parentMax < parentMin) {
                return false;
            }

            AmaxP = A->getUBounds()->at(pos);
            AminP = A->getLBounds()->at(pos);
            AmaxV = A->getUBounds()->at(vel);
            AminV = A->getLBounds()->at(vel);
        
            BmaxP = B->getUBounds()->at(pos);
            BminP = B->getLBounds()->at(pos);
            BmaxV = B->getUBounds()->at(vel);
            BminV = B->getLBounds()->at(vel);

                
            double newMinAcc, newMaxAcc;            
            double tmpAcc;
            std::vector<double> possibleAccs;
            
            /* calculate min acceleration first. set it to the highest
             * of the three values we could compute and the value
             * previously assigned. */
            possibleAccs.push_back(parentMin);
            
            tmpAcc = dt2*(BminP - AmaxP - AmaxV * dt);
            possibleAccs.push_back(tmpAcc);
                
            // jmyers: this item shaky - no one seems to understand it...
            tmpAcc = dt2*(AminP - BmaxP + BminV * dt);
            possibleAccs.push_back(tmpAcc);
                
            tmpAcc = (BminV - AmaxV) * dti;
            possibleAccs.push_back(tmpAcc);
                
            newMinAcc = *(std::max_element(possibleAccs.begin(),
                                           possibleAccs.end()));
            possibleAccs.clear();
                
            // now calculate new max acc.
            possibleAccs.push_back(parentMax);
                
            tmpAcc = dt2 * (BmaxP - AminP - AminV * dt);
            possibleAccs.push_back(tmpAcc);
                
            // jmyers: this item shaky - no one seems to understand it...
            tmpAcc = dt2 * (AmaxP - BminP + BmaxV * dt);
            possibleAccs.push_back(tmpAcc);
                
            tmpAcc = (BmaxV - AminV) * dti;
            possibleAccs.push_back(tmpAcc);
                
            newMaxAcc = *(std::min_element(possibleAccs.begin(),
                                           possibleAccs.end()));
            possibleAccs.clear();
                
            if (axis == 0) {
                minR = newMinAcc;
                maxR = newMaxAcc;
            }
            else {
                minD = newMinAcc;
                maxD = newMaxAcc;
            }
                
            // short-circuit if possible
            if (newMaxAcc < newMinAcc) {
                return false;
            }
        }
        
    }
    // we didn't short circuit so it must be valid. return true.
    return true;
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
            if (decDistance < searchConfig.trackAdditionThreshold) {

                double distance = angularDistanceRADec_deg(detRa, detDec, predRa, predDec);
                
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
 * Calculate track RMS. Trust that the best fit quadratic has been
 * calculated since the last addition of a detection and thus is up to
 * date.  Return true iff RMS < searchConfig.maxTrackRms.
 */ 
bool trackRmsIsSufficientlyLow(
    const std::vector<MopsDetection> &allDetections,
    const Track &newTrack, 
    const linkTrackletsConfig &searchConfig)
{

    double netSqError = 0.;
    std::set<uint> trackDets = newTrack.getComponentDetectionIndices();

    std::set<uint>::const_iterator detIter;
    for (detIter = trackDets.begin(); 
         detIter != trackDets.end(); 
         detIter++) {
        const MopsDetection thisDet = allDetections.at(*detIter);
        double predictedRa, predictedDec;
        newTrack.predictLocationAtTime(thisDet.getEpochMJD(), 
                                       predictedRa, 
                                       predictedDec);
        double dist = angularDistanceRADec_deg(thisDet.getRA(), 
                                               thisDet.getDec(), 
                                               predictedRa, 
                                               predictedDec);
        netSqError += dist*dist;
    }

    double rmsError = sqrt(netSqError / trackDets.size());

    return (rmsError < searchConfig.trackMaxRms);
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
    double start = std::clock();
    buildTracksAddToResultsVisits++;

    uint numCompatible = 0;

    //debugPrint(firstEndpoint, secondEndpoint, supportNodes,
    //allDetections, allTracklets);

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
            
            newTrack.calculateBestFitQuadratic(allDetections);

            if (endpointTrackletsAreCompatible(allDetections, 
                                               newTrack,
                                               searchConfig)) {
                
                numCompatible++;
                compatibleEndpointsFound++;
                
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
                    newTrack.calculateBestFitQuadratic(allDetections);
                    if (trackRmsIsSufficientlyLow(allDetections, 
                                                  newTrack, 
                                                  searchConfig)) {
                        results.insert(newTrack);
                    }
                }

                if ((RACE_TO_MAX_COMPATIBLE == true) && 
                    (compatibleEndpointsFound >= MAX_COMPATIBLE_TO_FIND)) {
                    debugPrintTimingInfo(results);
                    exit(0);
                }

            }
        }    
    }

    buildTracksAddToResultsTime += timeSince(start);

}






bool areAllLeaves(const std::vector<TreeNodeAndTime> &nodeArray) {
    double start = std::clock();
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
    areAllLeavesTime += timeSince(start);
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
    double start = std::clock();
    double width = 1;
    for (uint i = 0; i < 4; i++) {
        width *= node->getUBounds()->at(i) - node->getLBounds()->at(i);
    }
    nodeWidthTime += timeSince(start);
    return width;    
}





/* 
 * feb 17, 2011: update acc bounds using formulas reverse-engineered
 * from Kubica.
 */
bool updateAccBoundsReturnValidity(const TreeNodeAndTime &firstEndpoint, 
                                   const TreeNodeAndTime &secondEndpoint,
                                   double &accMinRa, double &accMaxRa, 
                                   double &accMinDec, double &accMaxDec)
{
    unsigned int axis;
    double parentMin, parentMax;

    double dt = secondEndpoint.myTime.getMJD() - 
        firstEndpoint.myTime.getMJD();
    double dt2 = 2./(dt*dt);
    double dti = 1./(dt);

    TrackletTreeNode* node1 = firstEndpoint.myTree;
    TrackletTreeNode* node2 = secondEndpoint.myTree;

    for (axis = 0; axis < 2; axis++) {
        double node1maxV, node1minV, node1maxP, node1minP;
        double node2maxV, node2minV, node2maxP, node2minP;
        unsigned int pos, vel;
        if (axis == 0) {
            pos=POINT_RA;
            vel=POINT_RA_VELOCITY;
            parentMin = accMinRa;
            parentMax = accMaxRa;
        }
        else {
            pos=POINT_DEC;
            vel=POINT_DEC_VELOCITY;
            parentMin = accMinDec;
            parentMax = accMaxDec;
        }

        // short-circuit right away if possible.
        if (parentMax < parentMin) return false;

        node1maxP = node1->getUBounds()->at(pos);
        node1minP = node1->getLBounds()->at(pos);
        node1maxV = node1->getUBounds()->at(vel);
        node1minV = node1->getLBounds()->at(vel);
        
        node2maxP = node2->getUBounds()->at(pos);
        node2minP = node2->getLBounds()->at(pos);
        node2maxV = node2->getUBounds()->at(vel);
        node2minV = node2->getLBounds()->at(vel);

        
        double newMinAcc, newMaxAcc;
            
        double tmpAcc;
        std::vector<double> possibleAccs;
        /* calculate min acceleration first. set it to the highest of
         * the three values we could compute and the value previously
         * assigned. */
        possibleAccs.push_back(parentMin);
            
        tmpAcc = (node2minV - node1maxV) * dti;
        possibleAccs.push_back(tmpAcc);
            
        tmpAcc = dt2 * (node2minP - node1maxP - node1maxV * dt);
        possibleAccs.push_back(tmpAcc);
            
        // jmyers: this item shaky - no one seems to understand it...
        tmpAcc = dt2 * (node1minP - node2maxP + node2minV * dt);
        possibleAccs.push_back(tmpAcc);
            
        newMinAcc = *(std::max_element(possibleAccs.begin(), 
                                       possibleAccs.end()));
        possibleAccs.clear();
            
        // now calculate max acceleration and take the least of the
        // possible values.
        possibleAccs.push_back(parentMax);
            
        tmpAcc = (node2maxV - node1minV) * dti;
        possibleAccs.push_back(tmpAcc);
            
        tmpAcc = dt2 * (node2maxP - node1minP - node1minV * dt);
        possibleAccs.push_back(tmpAcc);
            
        // jmyers: this item shaky - no one seems to understand it...
        tmpAcc = dt2 * (node1maxP - node2minP + node2maxV * dt);
        possibleAccs.push_back(tmpAcc);
            
        newMaxAcc = *(std::min_element(possibleAccs.begin(), 
                                       possibleAccs.end()));
            
        possibleAccs.clear();
            
        if (axis == 0) {
            accMinRa = newMinAcc;
            accMaxRa = newMaxAcc;
        }
        else {
            accMinDec = newMinAcc;
            accMaxDec = newMaxAcc;
        }
            
        // short-circuit if possible
        if (newMaxAcc < newMinAcc) return false;

    }
    // we know maxAcc > minAcc because we didn't short-circuit above.
    return true;
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
    //std::cout << " Starting doLinkingRecurse, acc bounds are:\n" 
    //          << accMinRa << " , " << accMaxRa << " ; " 
    //          << accMinDec << " , " << accMaxDec << ";\n";
    double start = std::clock();
    firstEndpoint.myTree->addVisit();

    doLinkingRecurseVisits++;

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
        }        

        if (iterationsTillSplit <= 0) {
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
                        //first endpoint.\n";
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
                        //left child of first endpoint.\n";
                    }
                    
                    if (firstEndpoint.myTree->hasRightChild())
                    {
                        TreeNodeAndTime newTAT(
                            firstEndpoint.myTree->getRightChild(), 
                            firstEndpoint.myTime);
                        doLinkingRecurseTime += timeSince(start);
                        //std::cout << "recursing on right child of
                        //first endpoint..\n";
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
                        //right child of first endpoint.\n";
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
                        doLinkingRecurseTime += timeSince(start);
                        //std::cout << "Recursing on left child of
                        //second endpoint.\n";
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
                        //left child of second endpoint.\n";
                    }
                    
                    if (secondEndpoint.myTree->hasRightChild())
                    {
                        TreeNodeAndTime newTAT(
                            secondEndpoint.myTree->getRightChild(), 
                            secondEndpoint.myTime);
                        doLinkingRecurseTime += timeSince(start);
                        //std::cout << "Recursing on right child of
                        //second endpoint.\n";
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
                        //right child of second endpoint.\n";
                        
                    }
                }
            }                        
        }
    }
}






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

    initDebugTimingInfo();

    unsigned int numImages = trackletTimeToTreeMap.size();

    std::map<ImageTime, TrackletTree >::const_iterator firstEndpointIter;

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
                
                        std::vector<TreeNodeAndTime > supportPoints;
                        std::map<ImageTime, TrackletTree >::const_iterator 
                            supportPointIter;

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
                
                        if (DEBUG) {std::cout << "\n";}

                        TreeNodeAndTime firstEndpoint(firstEndpointIter->second.getRootNode(), 
                                                      firstEndpointIter->first);
                        TreeNodeAndTime secondEndpoint(secondEndpointIter->second.getRootNode(),
                                                       secondEndpointIter->first);
                
                        //call the recursive linker with the endpoint
                        //nodes and support point nodes.

                        // use a cache to avoid doing math in
                        // isCompatible!  we use a new cache for each
                        // pair of endpoints, because isCompatible
                        // projects the location of the endpoint
                        // regions forward/backwards in time, so there
                        // will be 0 reuse between pairs of endpoint
                        // trees.

                        double iterationTime = std::clock();
                        if (searchConfig.myVerbosity.printStatus) {
                            std::cout << "Looking for tracks between images at times " 
                                      << std::setprecision(12) 
                                      << firstEndpointIter->first.getMJD() 
                                      << " (image " << firstEndpointIter->first.getImageId() 
                                      << " / " << numImages << ")"
                                      << " and " 
                                      << std::setprecision(12)
                                      << secondEndpointIter->first.getMJD()
                                      << " (image " 
                                      << secondEndpointIter->first.getImageId() 
                                      << " / " << numImages << ") "
                                      << " (with " 
                                      << supportPoints.size() << " support images).\n";
                            struct tm * timeinfo;
                            time_t rawtime;
                            time ( &rawtime );
                            timeinfo = localtime ( &rawtime );                    
                            
                            std::cout << " current wall-clock time is " 
                                      << asctime (timeinfo);

                        }
                        imagePairs += 1;
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
                                         results, 
                                         ITERATIONS_PER_SPLIT);

                        if (searchConfig.myVerbosity.printStatus) {
                            time_t rawtime;
                            
                            struct tm * timeinfo;
                            std::cout << "That iteration took " 
                                      << timeSince(iterationTime) << " seconds. "
                                      << std::endl;
                            std::cout << " so far, we have found " << 
                                results.size() << " tracks.\n\n";
                            time ( &rawtime );
                            timeinfo = localtime ( &rawtime );                    
                        }
                    }
                }
            }
        }
    }
    if (searchConfig.myVerbosity.printVisitCounts) {
        std::cout << "Found " << imagePairs << 
            " valid start/end image pairs.\n";
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
        std::cout << "Recentering all detections on (180, 0).\n";
    }
    recenterDetections(allDetections, searchConfig);
    
    if (searchConfig.myVerbosity.printStatus) {
        std::cout << "Setting tracklet velocities.\n";
    }
    setTrackletVelocities(allDetections, queryTracklets);

    if (searchConfig.myVerbosity.printStatus) {
        std::cout << "Sorting tracklets by image time and creating trees.\n";
    }
    std::map<ImageTime, TrackletTree > trackletTimeToTreeMap;    
    makeTrackletTimeToTreeMap(allDetections, 
                              queryTracklets, 
                              trackletTimeToTreeMap, 
                              searchConfig);
    if (searchConfig.myVerbosity.printStatus) {
        std::cout << "Doing the linking.\n";
    }

    clock_t linkingStart = std::clock();

    doLinking(allDetections, 
              queryTracklets, 
              searchConfig, 
              trackletTimeToTreeMap, 
              *toRet);
    if (searchConfig.myVerbosity.printStatus) {
        std::cout << "Finished linking.\n";
    }
    if (searchConfig.myVerbosity.printVisitCounts) {
        std::cout << "Made " << doLinkingRecurseVisits 
                  << " calls to doLinkingRecurse.\n";
        std::cout << "Made " << buildTracksAddToResultsVisits << 
            " calls to buildTracksAddToResults.\n";
        double linkingTime = timeSince(linkingStart);
        std::cout << "Linking took " << linkingTime << " seconds.\n";
        std::cout << " " << buildTracksAddToResultsTime 
                  << " sec spent on terminal tracklet processing and " 
                  << linkingTime - buildTracksAddToResultsTime
                  << " sec on other processing.\n";
        
    }
    
    return toRet;
}





}} //close lsst::mops

// -*- LSST-C++ -*-
/* jonathan myers */

// time headers needed for benchmarking performance
#include <ctime>
#include <iomanip>
#include <map>
#include <gsl/gsl_multifit.h>
#include <time.h>
#include <algorithm>

#include "lsst/mops/daymops/linkTracklets/linkTracklets.h"
#include "lsst/mops/rmsLineFit.h"
#include "lsst/mops/Exceptions.h"
#include "lsst/mops/KDTree.h"

#include "lsst/mops/daymops/linkTracklets/lruCache.h"



#define CACHE_SIZE 100000
#define USE_CACHE true

/* taking a queue from Kubica, it's only once per ITERATIONS_PER_SPLIT calls to
 * doLinkingRecurse that we actually split the (non-leaf) support nodes. The
 * idea is to avoid redundantly calculating whether a set of support nodes is
 * compatible.
 */
#define ITERATIONS_PER_SPLIT 0

/* 
the following flags, if set to 'true, will enable some debugging checks which
use brute-force searching and ground truth data to alert the user if an object
is going to be missed.  (not foolproof, but pretty good).  Will make searching
INSANELY slow.
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



namespace lsst {
    namespace mops {



double getAllDetectionsForTrackletTime, getFirstDetectionForTrackletTime, 
    modifyWithAccelerationTime, positionAndVelocityRangesOverlapAfterAccelerationTime, 
    areCompatibleTime, getBestFitVelocityAndAccelerationTime,
    getBestFitVelocityAndAccelerationForTrackletsTime, addBestCompatibleTrackletsAndDetectionsToTrackTime, 
    trackMeetsRequirementsTime, endpointTrackletsAreCompatibleTime, buildTracksAddToResultsTime, 
    areAllLeavesTime, nodeWidthTime, addAllDetectedObjectsToSetTime, doLinkingRecurseTime;

int doLinkingRecurseVisits, buildTracksAddToResultsVisits, compatibleEndpointsFound;
int rejectedOnVelocity, rejectedOnPosition, wereCompatible, rejectedOnLackOfSupport;
int cacheHits, cacheMisses;


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



class ImageTime {
public:
    ImageTime() {
        MJD = -1; imgId = 0;
    }
    ImageTime(double newMJD, unsigned int newImageId) {
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
    const unsigned int getImageId() const {
        return imgId;
    }
    void setMJD(double m) {
        MJD = m;
    }
    void setImageId(unsigned int i) {
        imgId = i;
    }
    ImageTime & operator=(const ImageTime &rhs) {
        MJD = rhs.getMJD();
        imgId = rhs.getImageId();
        return *this;
    }
    // NB: since IDs are assigned in order of image time, would it be faster if
    // we did < based on imgId (int rather than double?)
    bool operator< (const ImageTime &other) const {
        return MJD < other.getMJD();
    }

private:
    double MJD;
    unsigned int imgId;
};


/*
 * this class is only used as a key to the cache of TreeNode, img time T -> node
 * bounding box at time T. 
 *
 * note that KDTreeNodes have unique IDs *WITHIN THEIR TREE* but every tree will
 * contain a node with ID 1.
 *
 */
class GlobalNodeIdAndProjectionId {
public:
    GlobalNodeIdAndProjectionId(unsigned int _treeNodeId, 
                                unsigned int _treeImgTime, 
                                unsigned int _projectionTime) {
        treeNodeId = _treeNodeId;
        imgId = _treeImgTime;
        projectionId = _projectionTime;
    }
    GlobalNodeIdAndProjectionId() {
        treeNodeId = 0;
        imgId = 0;
    }

    void setImageId(unsigned int newId) {
        imgId = newId;
    }

    void setTreeNodeId(unsigned int newId) {
        treeNodeId = newId;
    }
    void setProjectionId(unsigned int newId) {
        projectionId = newId;
    }

    unsigned int getImageId() const {
        return imgId;
    }
    unsigned int getTreeNodeId() const {
        return treeNodeId;
    }
    unsigned int getProjectionId() const {
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
    unsigned int imgId;
    unsigned int treeNodeId;
    unsigned int projectionId;
};




// a quick macro to remove a lot of ugly typing
#define LTCache LRUCache <GlobalNodeIdAndProjectionId, std::vector <std::vector <double> > >




class TreeNodeAndTime {
public:
    TreeNodeAndTime(KDTreeNode<unsigned int> * tree, ImageTime i) {
        myTree = tree;
        myTime = i;
    }
    KDTreeNode <unsigned int> * myTree;
    ImageTime myTime;

};





std::set<unsigned int> setUnion(const std::set<unsigned int> &s1, 
                                const std::set<unsigned int> &s2) 
{
    std::set<unsigned int> toRet;
    std::set<unsigned int>::const_iterator iter;
    for (iter = s1.begin(); iter != s1.end(); iter++) {
        toRet.insert(*iter);
    }
    for (iter = s2.begin(); iter != s2.end(); iter++) {
        toRet.insert(*iter);
    }
    return toRet;
}






std::set<unsigned int> allDetsInTreeNode(KDTreeNode<unsigned int> &t,
                                            const std::vector<MopsDetection>&allDets,
                                            const std::vector<Tracklet>&allTracklets) 
{
    std::set<unsigned int> toRet;
    std::vector<PointAndValue<unsigned int> >::const_iterator tIter;

    if(! t.isLeaf() ) {
        if (t.hasLeftChild()) {
            std::set<unsigned int> childDets = allDetsInTreeNode(*t.getLeftChild(), allDets, allTracklets);
            toRet = setUnion(toRet, childDets);
        }
        if (t.hasRightChild()) {
            std::set<unsigned int> childDets = allDetsInTreeNode(*t.getRightChild(), allDets, allTracklets);
            toRet = setUnion(toRet, childDets);
        }
    }
    else 
    {
        // t.getMyData() holds points and values. the "value" part is index into allTracklets.
        for (tIter = t.getMyData()->begin();
             tIter != t.getMyData()->end();
             tIter++) {
            
            std::set<unsigned int>::const_iterator detIter;
            for (detIter = allTracklets.at(tIter->getValue()).indices.begin();
                 detIter != allTracklets.at(tIter->getValue()).indices.end();
                 detIter++) {
                
                toRet.insert(allDets.at(*detIter).getID());
                
            }
        }
    }

    return toRet;

}





void printSet(const std::set<unsigned int> s, std::string delimiter) 
{
    std::set<unsigned int>::const_iterator setIter;
    for (setIter = s.begin();
         setIter != s.end();
         setIter++) {
        std::cout << *setIter << delimiter;
    }

}








void debugPrint(const TreeNodeAndTime &firstEndpoint, const TreeNodeAndTime &secondEndpoint, 
                std::vector<TreeNodeAndTime> &supportNodes, 
                const std::vector<MopsDetection> &allDetections,
                const std::vector<Tracklet> &allTracklets) 
{
    std::set<unsigned int> leftEndpointDetIds = allDetsInTreeNode(*(firstEndpoint.myTree), 
                                                                  allDetections, allTracklets);
    std::set<unsigned int> rightEndpointDetIds = allDetsInTreeNode(*(secondEndpoint.myTree),
                                                                   allDetections, allTracklets);

    std::cout << "in doLinkingRecurse2,       first endpoint contains     " ;
    printSet(leftEndpointDetIds, " ");
    std::cout << '\n';
    std::cout << "                            second endpoint contains    " ;
    printSet(rightEndpointDetIds, " ");
    std::cout << '\n';

    for(unsigned int i = 0; i < supportNodes.size(); i++) {
        std::cout << " support node " << i << " contains                      ";
        std::set <unsigned int> supDetIds = allDetsInTreeNode(*(supportNodes.at(i).myTree), 
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




void showNumVisits(KDTreeNode<unsigned int> *tree, unsigned int &totalNodes, unsigned int &totalVisits) 
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




// the final parameter is modified; it will hold Detections associated with the 
// tracklet t.
void getAllDetectionsForTracklet(const std::vector<MopsDetection> & allDetections,
                              const Tracklet &t,
                              std::vector<MopsDetection> &detectionsForTracklet) 
{
    double start = std::clock();
    
    detectionsForTracklet.clear();
    std::set<unsigned int>::const_iterator trackletIndexIter;

    for (trackletIndexIter = t.indices.begin();
         trackletIndexIter != t.indices.end();
         trackletIndexIter++) {
        detectionsForTracklet.push_back(allDetections.at(*trackletIndexIter));
    }

    getAllDetectionsForTrackletTime += getTimeElapsed(start);
}






// given a tracklet, return the *temporally* earliest detection of this tracklet
// (the one with minimum MJD). 
MopsDetection getFirstDetectionForTracklet(const std::vector<MopsDetection> &allDetections,
                                           const Tracklet &t) 
{
    double start = std::clock();
    if (t.indices.size() < 1) {
        LSST_EXCEPT(BadParameterException,
                    "linkTracklets::getFirstDetectionForTracklet called with empty tracklet.");
    }

    MopsDetection toRet;
    bool foundOne = false;
    std::set<unsigned int>::const_iterator indexIter;
    for (indexIter = t.indices.begin(); indexIter != t.indices.end(); indexIter++) {
        MopsDetection curDet = allDetections.at(*indexIter);
        if ((!foundOne) || 
            (toRet.getEpochMJD() > curDet.getEpochMJD())) {
            toRet = curDet;
            foundOne = true;
        }
    }

    getFirstDetectionForTrackletTime += getTimeElapsed(start);
    return toRet;

}










void makeTrackletTimeToTreeMap(const std::vector<MopsDetection> &allDetections,
                               const std::vector<Tracklet> &queryTracklets,
                               std::map<ImageTime, KDTree <unsigned int> > &newMap,
                               linkTrackletsConfig myConf)
{
    bool printDebug = false;
    if (printDebug) {
        std::cout << "Sorting tracklets by \"root\" time. We have  " << allDetections.size() 
                  << " detections and " << queryTracklets.size() << " trakclets.\n";

    }

    newMap.clear();
    //sort all tracklets by their first image time; make PointAndValues from
    //these so we can build a tree.
    // allTrackletPAVsMap will map from image time -> [ all tracklets starting at that image time. ]
    std::map<double, std::vector<PointAndValue<unsigned int> > > allTrackletPAVsMap;

    allTrackletPAVsMap.clear();

    for (unsigned int i = 0; i < queryTracklets.size(); i++) {
        MopsDetection firstDetection =  getFirstDetectionForTracklet(allDetections, queryTracklets.at(i));
        double firstDetectionTime = firstDetection.getEpochMJD();
        PointAndValue<unsigned int> trackletPAV;
        std::vector<double> trackletPoint;

        trackletPoint.push_back(firstDetection.getRA());
        trackletPoint.push_back(firstDetection.getDec());
        trackletPoint.push_back(queryTracklets.at(i).velocityRA);
        trackletPoint.push_back(queryTracklets.at(i).velocityDec);
        trackletPAV.setPoint(trackletPoint);        

        trackletPAV.setValue(i);

        allTrackletPAVsMap[firstDetectionTime].push_back(trackletPAV);
    }

    if (printDebug) 
        std::cout << " got " << allTrackletPAVsMap.size() << " image times.\n";
    
    // iterate over each time/pointAndValueVec pair and build a corresponding
    // time/KDTree pair.  note that we iterate over a Map which uses MJD as key;
    // Maps sort their data by their key, so we are iterating over all image
    // times in order.
    unsigned int curImageID = 0;

    std::map<double, std::vector<PointAndValue<unsigned int> > >::const_iterator PAVIter;
    
    for (PAVIter = allTrackletPAVsMap.begin(); PAVIter != allTrackletPAVsMap.end(); PAVIter++) {


        KDTree<unsigned int> curTree(PAVIter->second, 4, myConf.leafSize);
        newMap[ImageTime(PAVIter->first, curImageID)] = curTree;

        if (printDebug) {
            std::cout << " image time " << PAVIter->first << " (with ID = " << curImageID << ") had " 
                      << PAVIter->second.size() << " tracklets, generating a tree of size " 
                      << curTree.size() << std::endl;
            std::cout << "    with leaf node size = " << myConf.leafSize << ", and an average leaf size of about " 
                      << PAVIter->second.size() * 2. / curTree.size() << std::endl;
        }

        curImageID++;
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
    double newPosition = position + velocity*time + .5*acceleration*(time*time);
    double newVelocity = velocity + acceleration*time;
    position = newPosition;
    velocity = newVelocity;
    modifyWithAccelerationTime += timeSince(start);
}










inline bool positionAndVelocityRangesOverlap(double firstPositionMin, double firstPositionMax, 
                                             double firstVelocityMin, double firstVelocityMax,
                                             double secondPositionMin, double secondPositionMax,
                                             double secondVelocityMin, double secondVelocityMax)
{

    double start = std::clock();

    bool velocityCompatible = regionsOverlap1D_unsafe(firstVelocityMin, firstVelocityMax,
                                                                      secondVelocityMin, secondVelocityMax);
    if (!velocityCompatible) {
        rejectedOnVelocity++;
        positionAndVelocityRangesOverlapAfterAccelerationTime += timeSince(start);
        return false;
    }

    firstPositionMin = convertToStandardDegrees(firstPositionMin);
    firstPositionMax = convertToStandardDegrees(firstPositionMax);

    secondPositionMin = convertToStandardDegrees(secondPositionMin);
    secondPositionMax = convertToStandardDegrees(secondPositionMax);

    bool positionCompatible = angularRegionsOverlapSafe(firstPositionMin, firstPositionMax,
                                                                        secondPositionMin, secondPositionMax);
    if (!positionCompatible) {
        rejectedOnPosition++;
        positionAndVelocityRangesOverlapAfterAccelerationTime += timeSince(start);
        return false;
    }
    
    wereCompatible++;
    positionAndVelocityRangesOverlapAfterAccelerationTime += timeSince(start);
    return true;
}



/*
 * to be used by areCompatible and probably nothing else.  extend the
 * position/velocity region BACKWARD in time; we expect deltaTime to be
 * NEGATIVE.
 *
 * This answers the question: what range contains all objects (accelerating no
 * faster than accel) that COULD reach this region of position/velocity space?
 */
void extendRangeBackward(double &p0Min,  double &p0Max,  double &vMin,
                         double &vMax,   double accel,  double deltaTime)
{
    if (deltaTime > 0) {
        throw LSST_EXCEPT(ProgrammerErrorException, 
                          "extendRangeBackward: Expected deltaTime to be negative!");
    }
    double newVMax = vMax + accel*fabs(deltaTime);
    double newVMin = vMin - accel*fabs(deltaTime);
    double newP0Max = p0Max + vMin*deltaTime + accel*deltaTime*deltaTime;
    double newP0Min = p0Min + vMax*deltaTime - accel*deltaTime*deltaTime;

    p0Min = newP0Min;
    p0Max = newP0Max;
    vMax = newVMax;
    vMin = newVMin;
}



void extendBoundsToTime(const KDTreeNode<unsigned int> * treeNode,
                        double deltaTime, 
                        const linkTrackletsConfig &searchConfig,
                        
                        const GlobalNodeIdAndProjectionId &lookupKey,
                        LTCache &rangeCache,
                        
                        double &raPositionMax,
                        double &raPositionMin,
                        double &raVelocityMax,
                        double &raVelocityMin,
                        
                        double &decPositionMax,
                        double &decPositionMin,
                        double &decVelocityMax,
                        double &decVelocityMin)
{
    std::vector<std::vector<double> > cachedBounds;

    //std::cout << "do we use cache?" << USE_CACHE << " cache has size " << rangeCache.size() << std::endl;
    if ((USE_CACHE) && (rangeCache.find(lookupKey, cachedBounds))) {
        //std::cout << " cache hit!\n";
        // we don't have to recompute it -it's there already!
        cacheHits++;
        std::vector<double> uBounds = cachedBounds.at(0);
        std::vector<double> lBounds = cachedBounds.at(1);
        
        raPositionMax = uBounds.at(POINT_RA);
        raPositionMin = lBounds.at(POINT_RA);

        raVelocityMax = uBounds.at(POINT_RA_VELOCITY);
        raVelocityMin = lBounds.at(POINT_RA_VELOCITY);
        
        decPositionMax = uBounds.at(POINT_DEC);
        decPositionMin = lBounds.at(POINT_DEC);

        decVelocityMax = uBounds.at(POINT_DEC_VELOCITY);
        decVelocityMin = lBounds.at(POINT_DEC_VELOCITY);
        
        if (
            (raVelocityMax < treeNode->getUBounds()->at(POINT_RA_VELOCITY))
            ||
            (raVelocityMin > treeNode->getLBounds()->at(POINT_RA_VELOCITY))) {
            throw LSST_EXCEPT(ProgrammerErrorException, 
                              "Found cached bounds at another time more restrictive than bounds at current time!");
        }
        
    }
    else {
        
        //std::cout << " cache miss...\n";
        // we have not yet projected this bounding box to the given image time. 
        // do so, adding in error, and save the results.
        
        cacheMisses++;
        
        //get upper and lower bounds of node's RA position, velocity
        //modify the positional bounds using error thresholds provided by user
        /*
          TBD: we need a way to track the possible extensions to the velocity range which would
          be caused by observational error.
        */
        
        
        raPositionMax = treeNode->getUBounds()->at(POINT_RA) + searchConfig.detectionLocationErrorThresh;
        raVelocityMax = treeNode->getUBounds()->at(POINT_RA_VELOCITY) + searchConfig.velocityErrorThresh;

        raPositionMin = treeNode->getLBounds()->at(POINT_RA) - searchConfig.detectionLocationErrorThresh;
        raVelocityMin = treeNode->getLBounds()->at(POINT_RA_VELOCITY) - searchConfig.velocityErrorThresh;

        //get upper and lower bounds of node's Dec position, velocity
                
        decPositionMax = treeNode->getUBounds()->at(POINT_DEC) + searchConfig.detectionLocationErrorThresh;
        decVelocityMax = treeNode->getUBounds()->at(POINT_DEC_VELOCITY) + searchConfig.velocityErrorThresh;

        decPositionMin = treeNode->getLBounds()->at(POINT_DEC) - searchConfig.detectionLocationErrorThresh;
        decVelocityMin = treeNode->getLBounds()->at(POINT_DEC_VELOCITY) - searchConfig.velocityErrorThresh;

        if (deltaTime >= 0) {

            modifyWithAcceleration(raPositionMax, raVelocityMax, 
                                   searchConfig.maxRAAccel, deltaTime);
            
            modifyWithAcceleration(raPositionMin, raVelocityMin, 
                                   searchConfig.maxRAAccel * -1.0, deltaTime);

            modifyWithAcceleration(decPositionMax, decVelocityMax, 
                                   searchConfig.maxDecAccel, deltaTime);
            
            modifyWithAcceleration(decPositionMin, decVelocityMin, 
                                   searchConfig.maxDecAccel * -1.0, fabs(deltaTime));
            
            
        }
        else {
            // negative delta time - we are asking what region of RA/Dec/dRA/dDec space could
            // *reach* this point 

            extendRangeBackward(raPositionMin,  raPositionMax,  
                                raVelocityMin,  raVelocityMax,  
                                searchConfig.maxRAAccel, deltaTime);

            extendRangeBackward(decPositionMin, decPositionMax, 
                                decVelocityMin, decVelocityMax, searchConfig.maxDecAccel, deltaTime);
            
        }

        if (USE_CACHE) {
            // save these newly-computed values to cache!        
            std::vector<std::vector <double> > boundsForStorage(2);
            std::vector<double> uBounds(4);
            std::vector<double> lBounds(4);
            
            uBounds.at(POINT_RA) =  raPositionMax;
            uBounds.at(POINT_DEC) = decPositionMax;
            uBounds.at(POINT_RA_VELOCITY) = raVelocityMax;
            uBounds.at(POINT_DEC_VELOCITY)= decVelocityMax;
            
            lBounds.at(POINT_RA) = raPositionMin;
            lBounds.at(POINT_DEC)= decPositionMin;
            lBounds.at(POINT_RA_VELOCITY) = raVelocityMin;
            lBounds.at(POINT_DEC_VELOCITY) = decVelocityMin;
            
            boundsForStorage.at(0) = uBounds;
            boundsForStorage.at(1) = lBounds;
            
            rangeCache.insert(lookupKey, boundsForStorage);
        }
    }
    
}



/*
 * a fairly mindless rewrite of the testing code in Kubica's
 * test_and_add_support_final from linker.c.
 * 
 * I don't fully understand the math I don't think that should matter!
 * 
 * we require that first, second, and third are in increasing chronological
 * order.
 */
bool areMutuallyCompatible(const TreeNodeAndTime &firstNode,
                           const TreeNodeAndTime &secondNode,
                           const TreeNodeAndTime &thirdNode,
                           linkTrackletsConfig searchConfig)
{
    bool valid = true;
    
    double maxR = searchConfig.maxRAAccel;
    double minR = searchConfig.maxRAAccel * -1.;
    double maxD = searchConfig.maxDecAccel;
    double minD = searchConfig.maxDecAccel * -1;
    
    // /* Test the support tree for validity against each model tree. */
    // for(i=0;(i<M)&&(valid);i++) {

    for (unsigned int whichPair = 0; whichPair < 2; whichPair++) {
        const KDTreeNode<unsigned int> * A;
        const KDTreeNode<unsigned int> * B;
        double aTime;
        double bTime;
        //   /* Load the model tree and compute the operation order. */
        //   mdl_tr = tbt_ptr_array_ref(mdl_pts,i);
        //   if(tbt_time(sup_tr) < tbt_time(mdl_tr)) {
        //     A = sup_tr; B = mdl_tr;
        //   } else {
        //     B = sup_tr; A = mdl_tr;
        //   }

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
    
        //   dt  = tbt_lo_time(B) - tbt_hi_time(A);
        //   dt2 = 2.0 / (dt * dt);
        //   dti = 1.0/dt;

        double dt = bTime - aTime;
        double dt2 = 2.0 / (dt * dt);
        double dti = 1.0/dt;
        double acc;

        double pe = searchConfig.detectionLocationErrorThresh;
        double ve = searchConfig.velocityErrorThresh;


        //   /* Do the velocity+position/position tests. */
        //   if(valid) {
    //  acc = dt2*  ((tbt_hi_RA(B)                         -  tbt_lo_RA(A)                       ) - tbt_lo_vRA(A)                                 * dt);
        acc = dt2 * (((B->getUBounds()->at(POINT_RA) + pe) - (A->getLBounds()->at(POINT_RA) - pe)) - (A->getLBounds()->at(POINT_RA_VELOCITY) - ve) * dt);
    //  if(maxR > acc) { maxR = acc; }
        if (maxR > acc) {maxR = acc; }
    //  acc = dt2*((tbt_lo_RA(B)                         -  tbt_hi_RA(A))                        -  tbt_hi_vRA(A)                                * dt);
        acc = dt2*(((B->getLBounds()->at(POINT_RA) - pe) - (A->getUBounds()->at(POINT_RA) + pe)) - (A->getUBounds()->at(POINT_RA_VELOCITY) + ve) * dt);
        //if(minR < acc) { minR = acc; }
        if (minR < acc)  {minR = acc; }
    //  acc  = dt2*(( tbt_hi_DEC(B)                        -  tbt_lo_DEC(A))                        -  tbt_lo_vDEC(A)                                * dt)
        acc = dt2 *(((B->getUBounds()->at(POINT_DEC) + pe) - (A->getLBounds()->at(POINT_DEC) - pe)) - (A->getLBounds()->at(POINT_DEC_VELOCITY) - ve) * dt);
    //  if(maxD > acc) { maxD = acc; }
        if (maxD > acc) { maxD = acc; }
    //  acc = dt2*(( tbt_lo_DEC(B)                       -  tbt_hi_DEC(A))                        -   tbt_hi_vDEC(A)                               * dt);
        acc = dt2*(((B->getLBounds()->at(POINT_DEC) -pe) - (A->getUBounds()->at(POINT_DEC) + pe)) - (A->getUBounds()->at(POINT_DEC_VELOCITY) + ve) * dt);
    //  if (minD < acc) { minD = acc; }
        if (minD < acc) {  minD = acc; }
    //  valid = (minD <= maxD)&&(minR <= maxR);
        valid = (minD <= maxD) && (minR <= maxR);
    
        ///* Do the velocity+position/position tests. */
        //     if(valid) {
        if (valid) {
        //  acc = dt2*(tbt_hi_RA(A)                         - tbt_lo_RA(B)                         + tbt_hi_vRA(B)                                *dt);
            acc = dt2*((A->getUBounds()->at(POINT_RA) + pe) - (B->getLBounds()->at(POINT_RA) - pe) + (B->getUBounds()->at(POINT_RA_VELOCITY) + ve)*dt);
        //  if(maxR > acc) { maxR = acc; }
            if(maxR > acc) { maxR = acc; }
        //  acc = dt2*( tbt_lo_RA(A)                        -   tbt_hi_RA(B)                       +  tbt_lo_vRA(B)                                * dt);
            acc = dt2*((A->getLBounds()->at(POINT_RA) - pe) - (B->getUBounds()->at(POINT_RA) + pe) + (B->getLBounds()->at(POINT_RA_VELOCITY) - ve) * dt);
        //  if(minR < acc) { minR = acc; }
            if(minR < acc) { minR = acc; }
        //  acc = dt2*( tbt_hi_DEC(A)                        -  tbt_lo_DEC(B)                        + tbt_hi_vDEC(B)                                 * dt);
            acc = dt2*((A->getUBounds()->at(POINT_DEC) + pe) - (B->getLBounds()->at(POINT_DEC) - pe) + (B->getUBounds()->at(POINT_DEC_VELOCITY) + ve) * dt);
        //  if(maxD > acc) { maxD = acc; }
            if(maxD > acc) { maxD = acc; }
        //  acc = dt2*( tbt_lo_DEC(A)                  -        tbt_hi_DEC(B)                        +  tbt_lo_vDEC(B)                               * dt);
            acc = dt2*((A->getLBounds()->at(POINT_DEC) - pe) - (B->getUBounds()->at(POINT_DEC) + pe) + (B->getLBounds()->at(POINT_DEC_VELOCITY) - ve) * dt);
        //  if(minD < acc) { minD = acc; }
            if(minD < acc) { minD = acc; }
        //  valid = (minD <= maxD)&&(minR <= maxR);
            valid = (minD <= maxD)&&(minR <= maxR);
            //       /* Determine the accel bounds with both velocity bounds. */
            //       if(valid) {
            if (valid) {
                //      acc = (tbt_hi_vRA(B)                         -  tbt_lo_vRA(A))                    * dti;
                acc = ((B->getUBounds()->at(POINT_RA_VELOCITY) + ve) - (A->getLBounds()->at(POINT_RA_VELOCITY) - ve)) * dti;
                //      if(maxR > acc) { maxR = acc; }
                if(maxR > acc) { maxR = acc; }
                //      acc = (tbt_lo_vRA(B)                    - tbt_hi_vRA(A)                   )*dti;
                acc = ((B->getLBounds()->at(POINT_RA_VELOCITY) - ve) - (A->getUBounds()->at(POINT_RA_VELOCITY) + ve))*dti;
                //      if(minR < acc) { minR = acc; }
                if(minR < acc) { minR = acc; }
                //      acc = (tbt_hi_vDEC(B)                          - tbt_lo_vDEC(A))                    *dti;
                acc = ((B->getUBounds()->at(POINT_DEC_VELOCITY) + ve) - (A->getLBounds()->at(POINT_DEC_VELOCITY) - ve)) *dti;
                //      if(maxD > acc) { maxD = acc; }
                if(maxD > acc) { maxD = acc; }
                //      acc = (tbt_lo_vDEC(B)                    - tbt_hi_vDEC(A))                    *dti;
                acc = ((B->getLBounds()->at(POINT_DEC_VELOCITY) - ve) - (A->getUBounds()->at(POINT_DEC_VELOCITY) + ve)) *dti;
                //      if(minD < acc) { minD = acc; }
                if(minD < acc) { minD = acc; }
                //      valid = (minR <= maxR)&&(minD <= maxD);
                valid = (minR <= maxR)&&(minD <= maxD);
            }
        }
    }

    return valid;
}







/*
  return whether two KDTreeNodes holding

  [RA, Dec, RAvelocity, Decvelocity] -> tracklet ID

  are "compatible": that is, whether an object contained in the first node could
  reach the second node.  This is computed using the upper and lower bounds
  on the nodes' [RA, Dec, RAvelocity, Decvelocity] values.

  Note that we treat RA, Dec as unrelated and euclidean.

  all tracklets held in the node are assumed to start at the specified times.

  are compatible according to the rules parameters provided by searchConfig.

 */



//TBD: searches should really be extended with quadraticErrorThresh
bool areCompatible(const TreeNodeAndTime  &nodeA,
                   const TreeNodeAndTime  &nodeB,
                   linkTrackletsConfig searchConfig,
                   LTCache &rangeCache)
{

    double firstRAPositionMax,  firstRAPositionMin,  secondRAPositionMax,  secondRAPositionMin;
    double firstDecPositionMax, firstDecPositionMin, secondDecPositionMax, secondDecPositionMin;
    double firstRAVelocityMax,  firstRAVelocityMin,  secondRAVelocityMax,  secondRAVelocityMin;
    double firstDecVelocityMax, firstDecVelocityMin, secondDecVelocityMax, secondDecVelocityMin;

    double start = std::clock();

    KDTreeNode<unsigned int> * first;
    double firstTime;

    KDTreeNode<unsigned int> * second;
    double secondTime;

    first = nodeA.myTree;
    firstTime = nodeA.myTime.getMJD();
    unsigned int firstTimeId = nodeA.myTime.getImageId();

    second = nodeB.myTree;
    secondTime = nodeB.myTime.getMJD();
    unsigned int secondTimeId = nodeB.myTime.getImageId();


    bool RACompatible = false;
    bool DecCompatible = false;

    // we want to project out to the MIDPOINT time between these two 
    // regions. This will have a smaller overlap than quadratically
    // extending one all the way out to the second time.
    double deltaTime = (secondTime - firstTime) / 2.;

    GlobalNodeIdAndProjectionId lookupKeyFirst(first->getId(), firstTimeId, secondTimeId);
    GlobalNodeIdAndProjectionId lookupKeySecond(second->getId(), secondTimeId, firstTimeId);


        
    extendBoundsToTime(first, deltaTime, searchConfig,
                       lookupKeyFirst, rangeCache,
                       firstRAPositionMax,  firstRAPositionMin,
                       firstRAVelocityMax,  firstRAVelocityMin,
                       firstDecPositionMax, firstDecPositionMin,
                       firstDecVelocityMax, firstDecVelocityMin);

    extendBoundsToTime(second, -1 * deltaTime, searchConfig,
                       lookupKeySecond, rangeCache,
                       secondRAPositionMax,  secondRAPositionMin,
                       secondRAVelocityMax,  secondRAVelocityMin,
                       secondDecPositionMax, secondDecPositionMin,
                       secondDecVelocityMax, secondDecVelocityMin);
                       
        
    
    DecCompatible = positionAndVelocityRangesOverlap(firstDecPositionMin, firstDecPositionMax,
                                                     firstDecVelocityMin, firstDecVelocityMax,
                                                     secondDecPositionMin, secondDecPositionMax,
                                                     secondDecVelocityMin, secondDecVelocityMax);
    if (!DecCompatible) {
        areCompatibleTime += timeSince(start);
        return false;

    }

    RACompatible = positionAndVelocityRangesOverlap(firstRAPositionMin, firstRAPositionMax,
                                                    firstRAVelocityMin, firstRAVelocityMax,
                                                    secondRAPositionMin, secondRAPositionMax,
                                                    secondRAVelocityMin, secondRAVelocityMax);


    if ((RACompatible == false)) {
        areCompatibleTime += timeSince(start);
        return false;
    }
    
    areCompatibleTime += timeSince(start);
    return true;
}









/*
  this is a potentially dangerous way of doing things.  we really need to know
  the min and max possible velocity for a tracklet, not the 'best fit' velocity
  for the tracklet.  (we can implicitly calculate the min and max if we assume
  2-point tracklets; but with many-point tracklets some additional possible
  error creeps in).

  TBD: be a tad more careful. for 2-point tracklets this is easy, for longer
  tracklets it may be trickier.
 */
void setTrackletVelocities(const std::vector<MopsDetection> &allDetections,
                           std::vector<Tracklet> &queryTracklets)
{
    for (unsigned int i = 0; i < queryTracklets.size(); i++) {
        Tracklet *curTracklet = &queryTracklets.at(i);
        std::vector <MopsDetection> trackletDets;
        getAllDetectionsForTracklet(allDetections, *curTracklet, trackletDets);

        std::vector<double> RASlopeAndOffset;
        std::vector<double> DecSlopeAndOffset;
        leastSquaresSolveForRADecLinear(&trackletDets,
                                                    RASlopeAndOffset,
                                                    DecSlopeAndOffset);
        
        curTracklet->velocityRA = RASlopeAndOffset.at(0);
        curTracklet->velocityDec = DecSlopeAndOffset.at(0);
    }

}






void getBestFitVelocityAndAcceleration(std::vector<double> positions, const std::vector<double> & times,
                                       double & velocity, double &acceleration, double &position0)
{
    if (positions.size() < 1) {
        velocity = 0; 
        acceleration = 0;
        position0 = 0;
        return;
    }

    // try to get these guys along the same 180-degree stretch of degrees...
    double p0 = positions.at(0);
    for (unsigned int i = 1; i < positions.size(); i++) {
        while ( positions.at(i) - p0 > 180) {
            positions.at(i) -= 360;
        }
        while ( p0 - positions.at(i) > 180) {
            positions.at(i) += 360;
        }
    }
    
    double start = std::clock();
    if (positions.size() != times.size()) {
        throw LSST_EXCEPT(ProgrammerErrorException,
                          "getBestFitVelocityAndAcceleration: position and time vectors not same size!");
    }

    /* we're using GSL for this. this is roughly adapted from the GSL
     * documentation; see
     * http://www.gnu.org/software/gsl/manual/html_node/Fitting-Examples.html
     */ 

    gsl_vector * y = gsl_vector_alloc(positions.size());
    gsl_matrix * X = gsl_matrix_alloc(positions.size(), 3);
        
    for (unsigned int i = 0; i < positions.size(); i++) {
        gsl_matrix_set(X, i, 0, 1.0);
        gsl_matrix_set(X, i, 1, times.at(i) );
        gsl_matrix_set(X, i, 2, times.at(i)*times.at(i) );
        gsl_vector_set(y, i, positions.at(i));
    }
    gsl_vector * c = gsl_vector_alloc(3); // times*c = positions 
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (positions.size(), 3);
    gsl_matrix * covariance = gsl_matrix_alloc(3,3);
    double chiSquared = 0;    
    // TBD: check return values for error, etc
    gsl_multifit_linear(X, y, c, covariance, &chiSquared, work);
    position0    = gsl_vector_get(c,0);
    velocity     = gsl_vector_get(c,1);
    acceleration = 2 * gsl_vector_get(c,2); // we need to multiply by 2 due to the 1/2 part of 1/2 at^2 

    gsl_vector_free(y);
    gsl_matrix_free(X);
    gsl_vector_free(c);
    gsl_matrix_free(covariance);
    gsl_multifit_linear_free(work);
    getBestFitVelocityAndAccelerationTime += timeSince(start);
}








void getBestFitVelocityAndAccelerationForTracklets(const std::vector<MopsDetection> &allDetections,
                                                   const std::vector<Tracklet> &queryTracklets,
                                                   const unsigned int trackletID1,
                                                   const unsigned int trackletID2,
                                                   double & RAVelocity, double & RAAcceleration, 
                                                   double & RAPosition0,
                                                   double & DecVelocity, double & DecAcceleration, 
                                                   double & DecPosition0, 
                                                   double & time0) 
{
    double start = std::clock();
    //get the detection IDs from tracklet 1 and tracklet 2 into a std::set.

    std::set <unsigned int> allDetectionIDs;
    std::set<unsigned int>::const_iterator trackletDetIter;
    for (trackletDetIter = queryTracklets.at(trackletID1).indices.begin(); 
         trackletDetIter != queryTracklets.at(trackletID1).indices.end();
         trackletDetIter++) {
        allDetectionIDs.insert(*trackletDetIter);
    }
    for (trackletDetIter = queryTracklets.at(trackletID2).indices.begin(); 
         trackletDetIter != queryTracklets.at(trackletID2).indices.end();
         trackletDetIter++) {
        allDetectionIDs.insert(*trackletDetIter);
    }

    // use a helper function to get the velocities and accelerations in RA, Dec.
    std::vector<double> RAs;
    std::vector<double> Decs;
    std::vector<double> times;
    std::set<unsigned int>::const_iterator detIter;

    double firstTime = allDetections.at(*allDetectionIDs.begin()).getEpochMJD();
    time0 = firstTime;
    
    for (detIter = allDetectionIDs.begin(); detIter != allDetectionIDs.end(); detIter++) {
        const MopsDetection * curDetection = &(allDetections.at(*detIter));
        RAs.push_back(curDetection->getRA());
        Decs.push_back(curDetection->getDec());
        times.push_back(curDetection->getEpochMJD() - firstTime);
    }
    getBestFitVelocityAndAcceleration(RAs, times,  RAVelocity,  RAAcceleration,  RAPosition0 );
    getBestFitVelocityAndAcceleration(Decs, times, DecVelocity, DecAcceleration, DecPosition0);
    getBestFitVelocityAndAccelerationForTrackletsTime += timeSince(start);
}


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
 * addBestCompatibleTrackletsAndDetectionsToTrack: this function uses an RA/Dec
 * motion prediction and finds the most compatible points owned by tracklets in
 * the candidateTrackletsIDs vector.  The best-fit points are added in order of
 * fit quality, with at most one detection added per image time.  "Compatible"
 * here is defined by the parameters in searchConfig.
 *
 * Detection IDs and the IDs of the detections' parents are added to newTrack's
 * relevant fields.
 */
void addBestCompatibleTrackletsAndDetectionsToTrack(const std::vector<MopsDetection> &allDetections, 
                                                    const std::vector<Tracklet> &allTracklets, 
                                                    const std::vector<unsigned int> candidateTrackletIDs, 
                                                    double RAVelocity, double RAAcceleration, double RAPosition0,
                                                    double DecVelocity, double DecAcceleration, double DecPosition0,
                                                    double time0,
                                                    const linkTrackletsConfig &searchConfig,
                                                    Track &newTrack) 
{
    double start = std::clock();

    /* the scoreToIDsMap will hold detection MJDs as the key, map from detection
     * time to best-fitting candidate detection at that time
     */
    std::map<double, CandidateDetection > timeToCandidateMap;

    // find the best compatible detection at each unique image time
    std::vector<unsigned int>::const_iterator trackletIDIter;
    std::set<unsigned int>::const_iterator detectionIDIter;
    for (trackletIDIter = candidateTrackletIDs.begin();
         trackletIDIter != candidateTrackletIDs.end();
         trackletIDIter++) {
        const Tracklet * curTracklet = &allTracklets.at(*trackletIDIter);
        for (detectionIDIter =  curTracklet->indices.begin();
             detectionIDIter != curTracklet->indices.end();
             detectionIDIter++) {
            double detMJD = allDetections.at(*detectionIDIter).getEpochMJD();
            double detRA  = allDetections.at(*detectionIDIter).getRA();
            double detDec = allDetections.at(*detectionIDIter).getDec();
            double timeOffset = detMJD - time0;
            double predRA  = RAPosition0 + RAVelocity*timeOffset + .5 * RAAcceleration*timeOffset*timeOffset;
            double predDec = DecPosition0 + DecVelocity*timeOffset + .5 * DecAcceleration*timeOffset*timeOffset;

            double decDistance = fabs(detDec - predDec);
            // filter on dec distance, which is easy to compute, before trying the angular distance.
            if (decDistance < searchConfig.quadraticFitErrorThresh + searchConfig.detectionLocationErrorThresh) {
                double distance = angularDistanceRADec_deg(detRA, detDec, predRA, predDec);
            
                // if the detection is compatible, consider whether it's the best at the image time
                if (distance < searchConfig.quadraticFitErrorThresh + searchConfig.detectionLocationErrorThresh) {
                    std::map<double, CandidateDetection>::iterator candidateAtTime;
                    candidateAtTime = timeToCandidateMap.find(detMJD);
                    
                    // more crazy C++-talk for "if we have no candidate, or if this is a
                    // better candidate than the one we have, then add this as the new
                    // candidate at that time"
                    if ((candidateAtTime == timeToCandidateMap.end()) 
                        || (candidateAtTime->second.distance > distance)) {
                        CandidateDetection newCandidate(distance, *detectionIDIter, *trackletIDIter);
                        timeToCandidateMap[detMJD] = newCandidate;
                    }
                }
            }
        }
    }
    
    /* initialize a list of image times present in the track already. */
    std::set<double> trackMJDs;
    std::set<unsigned int>::const_iterator trackDetectionIndices;
    for (trackDetectionIndices =  newTrack.componentDetectionIndices.begin();
         trackDetectionIndices != newTrack.componentDetectionIndices.end();
         trackDetectionIndices++) {
        trackMJDs.insert(allDetections.at(*trackDetectionIndices).getEpochMJD());        
    }
    
    /* add detections (and their parent tracklets) in order of 'score' (distance
     * from best-fit line) without adding any detections from
     * already-represented image times
     */
    std::map<double, CandidateDetection>::iterator candidatesIter;

    for (candidatesIter = timeToCandidateMap.begin(); 
         candidatesIter != timeToCandidateMap.end(); 
         candidatesIter++) {

        if (trackMJDs.find(candidatesIter->first) == trackMJDs.end()) {
            /* add this detection and tracklet to the track and add the associated MJD to 
             * trackMJDs */
            
            trackMJDs.insert(candidatesIter->first);
            newTrack.componentDetectionIndices.insert(candidatesIter->second.detId);
            newTrack.componentTrackletIndices.insert(candidatesIter->second.parentTrackletId);
        }
    }
    addBestCompatibleTrackletsAndDetectionsToTrackTime += timeSince(start);
}







bool trackMeetsRequirements(const std::vector<MopsDetection> & allDetections, 
                            const Track &newTrack, 
                            double RAVelocity, double RAAcceleration, double RAPosition0,
                            double DecVelocity, double DecAcceleration, double DecPosition0,
                            double time0,
                            linkTrackletsConfig searchConfig)
{
    double start = std::clock();
    bool meetsReqs = true;

    if (newTrack.componentTrackletIndices.size() < searchConfig.minSupportTracklets + 2) {
        meetsReqs = false;
    }

    if (newTrack.componentDetectionIndices.size() < searchConfig.minDetectionsPerTrack) {
        meetsReqs = false;
    }

    if ((RAAcceleration > searchConfig.maxRAAccel) || (DecAcceleration > searchConfig.maxDecAccel)) {
        meetsReqs = false;
    }

    trackMeetsRequirementsTime += timeSince(start);
    return meetsReqs;
    
    
}







/*
 * checks that endpoint tracklets are compatible.  
 * 
 * this checks several things and returns true iff all are true:
 *
 * - all points are within error threshholds of the predicted location (as predicected by *Velocity, *Acceleration *Position0).
 * - the first and last detection are at least minTimeSeparation apart 
 * - the best-fit accelerations are within min/max bounds
 * 
 */
bool endpointTrackletsAreCompatible(const std::vector<MopsDetection> & allDetections, 
                                    const std::vector<Tracklet> &allTracklets,
                                    unsigned int trackletID1,
                                    unsigned int trackletID2,
                                    double RAVelocity, double RAAcceleration, double RAPosition0,
                                    double DecVelocity, double DecAcceleration, double DecPosition0,
                                    double time0,
                                    linkTrackletsConfig searchConfig)
{
    double start = std::clock();
    bool allOK = true;

    if (RAAcceleration > searchConfig.maxRAAccel) {
        allOK = false;
    }
    if (DecAcceleration > searchConfig.maxDecAccel) {
        allOK = false;
    }

    if (allOK == true) {
        
        std::set<unsigned int> allDetectionIndices; // the set of all detections held by both tracklets
        std::set<unsigned int>::const_iterator detIter;
        for (detIter = allTracklets.at(trackletID1).indices.begin();
             detIter != allTracklets.at(trackletID1).indices.end();
             detIter++) {
            allDetectionIndices.insert(*detIter);
        }
        for (detIter = allTracklets.at(trackletID2).indices.begin();
             detIter != allTracklets.at(trackletID2).indices.end();
             detIter++) {
            allDetectionIndices.insert(*detIter);
        }
        
        
        for (detIter = allDetectionIndices.begin();
             detIter != allDetectionIndices.end();
             detIter++) {
            double detMJD = allDetections.at(*detIter).getEpochMJD();
            double timeOffset = detMJD - time0;
            double RAPred = RAPosition0 + RAVelocity*timeOffset + .5 * RAAcceleration*timeOffset*timeOffset;
            double DecPred = DecPosition0 + DecVelocity*timeOffset + .5 * DecAcceleration*timeOffset*timeOffset;
            double observedRA = allDetections.at(*detIter).getRA();
            double observedDec = allDetections.at(*detIter).getDec();
            double distanceError = 
                angularDistanceRADec_deg(RAPred, DecPred, observedRA, observedDec);
            if (distanceError > searchConfig.quadraticFitErrorThresh + searchConfig.detectionLocationErrorThresh) {
                allOK = false;
            }
        }
        
        
        if (allOK == true) {
            //check that time separation is good
            double minMJD, maxMJD;
            detIter = allDetectionIndices.begin();
            minMJD = allDetections.at(*detIter).getEpochMJD();
            maxMJD = minMJD;
            for (detIter = allDetectionIndices.begin();
                 detIter != allDetectionIndices.end();
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
    }
    endpointTrackletsAreCompatibleTime += timeSince(start);
    return allOK;
}








/*
 * this is called when all endpoint nodes (i.e. model nodes) and support nodes
 * are leaves.  model nodes and support nodes are expected to be mutually compatible.
 */
void buildTracksAddToResults(const std::vector<MopsDetection> &allDetections,
                             const std::vector<Tracklet> &allTracklets,
                             linkTrackletsConfig searchConfig,
                             TreeNodeAndTime &firstEndpoint,
                             TreeNodeAndTime &secondEndpoint,
                             std::vector<TreeNodeAndTime> &supportNodes,
                             TrackSet & results)
{
    double start = std::clock();
    buildTracksAddToResultsVisits++;

    unsigned int numCompatible = 0;
    double addBestStartTime;
    double addBestTotalTime = 0;


    //debugPrint(firstEndpoint, secondEndpoint, supportNodes, allDetections, allTracklets);

    if ((firstEndpoint.myTree->isLeaf() == false) ||
        (secondEndpoint.myTree->isLeaf() == false)) {
        LSST_EXCEPT(ProgrammerErrorException, 
                    "buildTracksAddToResults got non-leaf nodes, must be a bug!");
    }
    for (unsigned int i = 0; i < supportNodes.size(); i++) {
        if (supportNodes.at(i).myTree->isLeaf() == false) {
            LSST_EXCEPT(ProgrammerErrorException, 
                        "buildTracksAddToResults got non-leaf nodes, must be a bug!");            
        }
    }

    std::vector<PointAndValue<unsigned int> >::const_iterator firstEndpointIter;
    std::vector<PointAndValue<unsigned int> >::const_iterator secondEndpointIter;
    std::vector<TreeNodeAndTime>::const_iterator supportNodeIter;
    std::vector<PointAndValue <unsigned int> >::const_iterator supportPointIter;
    for (firstEndpointIter = firstEndpoint.myTree->getMyData()->begin();
         firstEndpointIter != firstEndpoint.myTree->getMyData()->end();
         firstEndpointIter++) {

        for (secondEndpointIter = secondEndpoint.myTree->getMyData()->begin();
             secondEndpointIter != secondEndpoint.myTree->getMyData()->end();
             secondEndpointIter++) {
            
            std::set<unsigned int>::const_iterator notFound = 
                allTracklets.at(firstEndpointIter->getValue()).indices.end();
            std::set<unsigned int>::const_iterator notFound2 = 
                allTracklets.at(secondEndpointIter->getValue()).indices.end();

            /* figure out the rough quadratic track fitting the two endpoints.
             * if error is too large, quit. Otherwise, choose support points
             * from the support nodes, using best-fit first, and ignoring those
             * too far off the line.  If we get enough points, return a track.
            */

            double RAVelocity, DecVelocity, RAAcceleration, DecAcceleration;
            double RAPosition0, DecPosition0;
            double time0;
            getBestFitVelocityAndAccelerationForTracklets(allDetections,
                                                          allTracklets, 
                                                          firstEndpointIter->getValue(),
                                                          secondEndpointIter->getValue(),
                                                          RAVelocity,RAAcceleration,RAPosition0,
                                                          DecVelocity,DecAcceleration,DecPosition0,
                                                          time0);

            if (endpointTrackletsAreCompatible(allDetections, 
                                               allTracklets,
                                               firstEndpointIter->getValue(),
                                               secondEndpointIter->getValue(),
                                               RAVelocity,RAAcceleration,RAPosition0,
                                               DecVelocity,DecAcceleration,DecPosition0,
                                               time0,searchConfig)) {
                
                numCompatible++;
                compatibleEndpointsFound++;
                // create a new track with these endpoints
                Track newTrack;
                std::vector<unsigned int> candidateTrackletIDs;
                candidateTrackletIDs.push_back(firstEndpointIter->getValue() );
                candidateTrackletIDs.push_back(secondEndpointIter->getValue());
                addBestCompatibleTrackletsAndDetectionsToTrack(allDetections, allTracklets, candidateTrackletIDs, 
                                                               RAVelocity, RAAcceleration, RAPosition0,
                                                               DecVelocity, DecAcceleration, DecPosition0,
                                                               time0,searchConfig,
                                                               newTrack);
                
                candidateTrackletIDs.clear();
                
                /* add the best compatible support points */

                // put all support tracklet IDs in curSupportNodeData, then call
                // addBestCompatibleTrackletsAndDetectionsToTrack
                for (supportNodeIter = supportNodes.begin(); supportNodeIter != supportNodes.end();
                     supportNodeIter++) {
                    const std::vector<PointAndValue <unsigned int> > * curSupportNodeData;
                    if (!supportNodeIter->myTree->isLeaf()) {
                        throw LSST_EXCEPT(BadParameterException,
                                          std::string(__FUNCTION__) + 
                                          std::string(": received non-leaf node as support node."));
                    }
                    curSupportNodeData = supportNodeIter->myTree->getMyData(); 
                    for (supportPointIter  = curSupportNodeData->begin(); 
                         supportPointIter != curSupportNodeData->end();
                         supportPointIter++) {
                        candidateTrackletIDs.push_back(supportPointIter->getValue());
                    }
                }
                addBestStartTime = std::clock();
                addBestCompatibleTrackletsAndDetectionsToTrack(allDetections, allTracklets, candidateTrackletIDs,
                                                               RAVelocity, RAAcceleration, RAPosition0,
                                                               DecVelocity, DecAcceleration, DecPosition0, time0,
                                                               searchConfig,
                                                               newTrack);
                addBestTotalTime += timeSince(addBestStartTime);
                
                if (trackMeetsRequirements(allDetections, newTrack,  
                                           RAVelocity, RAAcceleration, RAPosition0, 
                                           DecVelocity, DecAcceleration, DecPosition0,
                                           time0,searchConfig)) {
                    results.insert(newTrack);
                }

                if ((RACE_TO_MAX_COMPATIBLE == true) && (compatibleEndpointsFound >= MAX_COMPATIBLE_TO_FIND)) {
                    debugPrintTimingInfo(results);
                    exit(0);
                }

            }
        }    
    }

    /*std::cout << "BuildTracksAddToResults found " << numCompatible << " compatible endpoint and searched " <<
        supportNodes.size() << " support nodes for support.\n";
    std::cout << " this took "  << timeSince(start) << " sec, with " << addBestTotalTime
    << " spent on choosing support detections." << std::endl;*/
    buildTracksAddToResultsTime += timeSince(start);

}






bool areAllLeaves(const std::vector<TreeNodeAndTime> &nodeArray) {
    double start = std::clock();
    bool allLeaves = true;
    std::vector<TreeNodeAndTime>::const_iterator treeIter;
    unsigned int count = 0;
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







bool supportTooWide(const TreeNodeAndTime& firstEndpoint, const TreeNodeAndTime& secondEndpoint,
                    const TreeNodeAndTime& supportNode) 
{
    /* odd-looking stuff with alpha based on test_and_add_support in
     * Kubica's linker.c.  the idea is to weight the expected size of a
     * node according to its proximity to either the first or last
     * endpoint . */
    double alpha = (supportNode.myTime.getMJD() - firstEndpoint.myTime.getMJD()) /
        (secondEndpoint.myTime.getMJD() - firstEndpoint.myTime.getMJD());
    
    // check the width of the support node in all 4 axes; compare with width of 
    for (unsigned int i = 0; i < 4; i++) {
        double nodeWidth = supportNode.myTree->getUBounds()->at(i) - supportNode.myTree->getLBounds()->at(i);
        double maxWidth = (1.0 - alpha) * 
            (firstEndpoint.myTree->getUBounds()->at(i) - firstEndpoint.myTree->getLBounds()->at(i)) + 
            alpha * (secondEndpoint.myTree->getUBounds()->at(i) - secondEndpoint.myTree->getLBounds()->at(i));
        
        // this constant 4 is taken from Kubica's linker.c test_and_add_support. 
        // Ask him to justify it, not me!
        if (4.0 * nodeWidth > maxWidth) {
            return true;
        }
    }

    return false;
}







void splitSupportRecursively(const TreeNodeAndTime& firstEndpoint, const TreeNodeAndTime& secondEndpoint, 
                             bool requireLeaves,
                             const TreeNodeAndTime &supportNode, linkTrackletsConfig searchConfig, 
                             LTCache &rangeCache,
                             std::vector<TreeNodeAndTime> &newSupportNodes)
{

    if ((firstEndpoint.myTime.getMJD() >= supportNode.myTime.getMJD()) || 
        (supportNode.myTime.getMJD() > secondEndpoint.myTime.getMJD())) {
        throw LSST_EXCEPT(BadParameterException, "splitSupportRecursively got impossibly-ordered endpoints/support");
    }

    if ((areCompatible(firstEndpoint, supportNode, searchConfig, rangeCache) && 
       areCompatible(secondEndpoint, supportNode, searchConfig, rangeCache))) {
        //if (areMutuallyCompatible(firstEndpoint, supportNode, secondEndpoint, searchConfig)) {

        if (supportNode.myTree->isLeaf()) {
            newSupportNodes.push_back(supportNode);
        }

        else if (requireLeaves) {
            if (supportNode.myTree->hasLeftChild()) {
                TreeNodeAndTime leftTat(supportNode.myTree->getLeftChild(), supportNode.myTime); 
                splitSupportRecursively(firstEndpoint, secondEndpoint, requireLeaves, 
                                        leftTat,
                                        searchConfig, rangeCache, newSupportNodes);
            }
            if (supportNode.myTree->hasRightChild()) {
                TreeNodeAndTime rightTat(supportNode.myTree->getRightChild(), supportNode.myTime);
                splitSupportRecursively(firstEndpoint, secondEndpoint, requireLeaves, 
                                        rightTat,
                                        searchConfig, rangeCache, newSupportNodes);
            }
        }
        
        else {
            // we don't require leaves in output, but check to see if we *should* split this node.
            // if not, add it to output. Otherwise, recurse on its children.

            bool tooWide = supportTooWide(firstEndpoint, secondEndpoint, supportNode);
            
            if (tooWide) {
                if (supportNode.myTree->hasLeftChild()) {
                    TreeNodeAndTime leftTat(supportNode.myTree->getLeftChild(), supportNode.myTime); 
                    splitSupportRecursively(firstEndpoint, secondEndpoint, requireLeaves, 
                                            leftTat,
                                            searchConfig, rangeCache, newSupportNodes);
                }
                if (supportNode.myTree->hasRightChild()) {
                    TreeNodeAndTime rightTat(supportNode.myTree->getRightChild(), supportNode.myTime);
                    splitSupportRecursively(firstEndpoint, secondEndpoint, requireLeaves, 
                                            rightTat,
                                            searchConfig, rangeCache, newSupportNodes);
                }
            }
            else {
                newSupportNodes.push_back(supportNode);
            }
        }
    }
}





void filterAndSplitSupport(const TreeNodeAndTime& firstEndpoint, const TreeNodeAndTime& secondEndpoint, 
                           const std::vector<TreeNodeAndTime> &supportNodes, linkTrackletsConfig searchConfig, 
                           LTCache &rangeCache,
                           std::vector<TreeNodeAndTime> &newSupportNodes) 
{

    // if the endpoints are leaves, require that we get all leaves in the support nodes.
    bool endpointsAreLeaves = firstEndpoint.myTree->isLeaf() && secondEndpoint.myTree->isLeaf();
    
    for (unsigned int i = 0; i < supportNodes.size(); i++) {
        splitSupportRecursively(firstEndpoint, secondEndpoint, 
                                endpointsAreLeaves, // 'true' here means we require leaves in output
                                supportNodes[i],
                                searchConfig, rangeCache, newSupportNodes);
    }
    
    
}







double nodeWidth(KDTreeNode<unsigned int> *node)
{
    double start = std::clock();
    double width = 1;
    for (unsigned int i = 0; i < 4; i++) {
        width *= node->getUBounds()->at(i) - node->getLBounds()->at(i);
    }
    nodeWidthTime += timeSince(start);
    return width;    
}






/*
 * this is, roughly, the algorithm presented in http://arxiv.org/abs/astro-ph/0703475v1:
 * 
 * Efficient intra- and inter-night linking of asteroid detections using kd-trees
 * 
 * The "proper" version of the algorithm (as presented in psuedocode om the
 * document) is a little different; basically, he has us recurring only on
 * endpoints, and splitting support nodes repeatedly at each call, until they
 * are no longer "too wide".  Unfortunately, I found that neither my intuition
 * nor Kubica's implementation made clear the definition of "too wide".
 *
 * this implementation is more like the one Kubica did in his code. at every
 * step, we check all support nodes for compatibility, splitting each one. we
 * then split one model node and recurse. 
 */
void doLinkingRecurse2(const std::vector<MopsDetection> &allDetections,
                       const std::vector<Tracklet> &allTracklets,
                       linkTrackletsConfig searchConfig,
                       TreeNodeAndTime &firstEndpoint,
                       TreeNodeAndTime &secondEndpoint,
                       std::vector<TreeNodeAndTime> &supportNodes,
                       TrackSet & results,
                       int iterationsTillSplit,
                       LTCache &rangeCache)
{
    double start = std::clock();
    doLinkingRecurseVisits++;
    firstEndpoint.myTree->addVisit();
    //std::cout << "entering doLinkingRecurse2." << std::endl;
    //debugPrint(firstEndpoint, secondEndpoint, supportNodes, allDetections, allTracklets);

    //debugPrintTimingInfo(results);
    //std::cout << "\ndoLinkingRecurse called with endpoint IDs " << firstEndpoint.myTree->getId() << ",  " <<
    //    secondEndpoint.myTree->getId() << std::endl;
    //std::cout << " first endpoint " << (firstEndpoint.myTree->isLeaf() ? std::string("IS") : std::string("IS NOT")) << " a leaf; second endpoint " << 
    //    (secondEndpoint.myTree->isLeaf() ? std::string("IS") : std::string("IS NOT")) << " a leaf." << "\n";
    //std::cout << "   also we have " << supportNodes.size() << " support nodes." << std:: endl;
    
    if (areCompatible(firstEndpoint, secondEndpoint, searchConfig, rangeCache) == false)
    {
        // poor choice of model nodes (endpoint nodes)! give up.
        doLinkingRecurseTime += timeSince(start);
        return;
    }
    else 
    {

        std::set<double> uniqueSupportMJDs;
        std::vector<TreeNodeAndTime> newSupportNodes;
        std::vector<TreeNodeAndTime>::iterator supportNodeIter;
        
        /* look through untested support nodes, find the ones that are compatible
         * with the model nodes, add their children to newSupportNodes */

        if ((iterationsTillSplit <= 0) || 
            (firstEndpoint.myTree->isLeaf() && secondEndpoint.myTree->isLeaf())) {
            
            double filterStartTime = std::clock();
            //std::cout << "calling filterAndSplitSupport with " << supportNodes.size() << " support nodes.\n";
            filterAndSplitSupport(firstEndpoint, secondEndpoint, 
                                  supportNodes, searchConfig, rangeCache,
                                  newSupportNodes);
            //std::cout << " returned after " << timeSince(filterStartTime) << " and we now have " << 
            //    newSupportNodes.size() << " support nodes.\n";
        }        

        if (iterationsTillSplit <= 0) {
            iterationsTillSplit = ITERATIONS_PER_SPLIT;
        }
        else{
            // we still need to get newSupportNodes set up. just use the old ones.
            newSupportNodes = supportNodes;
        }
        /*std::cout << "Filtering support nodes took " << timeSince(start) << " sec.\n";
        std::cout << "After filtering support nodes, we now have " << newSupportNodes.size() 
        <<  " support nodes." << std::endl;*/
        
        std::vector<TreeNodeAndTime>::const_iterator supportIter;
        for (supportIter = newSupportNodes.begin(); supportIter != newSupportNodes.end(); supportIter++) {
            uniqueSupportMJDs.insert(supportIter->myTime.getMJD());
        }

        if (uniqueSupportMJDs.size() < searchConfig.minSupportTracklets) {
            // we can't possibly have enough distinct support tracklets between endpoints
            rejectedOnLackOfSupport++;
            doLinkingRecurseTime += timeSince(start);
            // std::cout << "Exiting doLinkingRecurse2 due to lack of support.\n";
            return; 
        }
        else 
        {
            // we have enough model nodes, and enough support nodes.
            // if they are all leaves, then start building tracks.
            // if they are not leaves, split one of them and recurse.

            if (firstEndpoint.myTree->isLeaf() && secondEndpoint.myTree->isLeaf() &&
                areAllLeaves(newSupportNodes)) {
                
                // TBD: actually, we need to filter newSupportNodes one more time here, since 
                // we checked for compatibility, then split off the children. of course, buildTracksAddToResults won't 
                // be affected with regards to correctness, just performance. So it may not even be wise to add a
                // mostly-needless check.
                //std::cout << "Calling buildTracksAddToResults\n";
                //double buildTracksStart = std::clock();
                //unsigned int numResults = results.size();
                buildTracksAddToResults(allDetections, allTracklets, searchConfig,
                                        firstEndpoint, secondEndpoint, newSupportNodes,
                                        results);
                //std::cout << "Brute-force track building took " << timeSince(buildTracksStart) << " sec.";
                //std::cout << "Found " << results.size() - numResults << " new tracks. Exiting." << std::endl;
            }
            else  {
               
                iterationsTillSplit -= 1;

                // find the "widest" node, where width is just the product
                // of RA range, Dec range, RA velocity range, Dec velocity range.
                // we will split that node and recurse.
                
                double firstEndpointWidth = nodeWidth(firstEndpoint.myTree);
                double secondEndpointWidth = nodeWidth(secondEndpoint.myTree);
                
                // don't consider splitting endpoint nodes which are actually leaves! give them negative width to hack selection process.
                if (firstEndpoint.myTree->isLeaf()) {
                    firstEndpointWidth = -1;
                }
                if (secondEndpoint.myTree->isLeaf()) {
                    secondEndpointWidth = -1;
                }
 
                // choose the widest model node, split it and recurse!                                
                if ( (areEqual(firstEndpointWidth, -1)) &&
                     (areEqual(secondEndpointWidth, -1)) ) {
                    // in this case, our endpoints are leaves, but not all our support nodes are.
                    // just call this function again until they *are* all leaves.
                    iterationsTillSplit = 0;
                    doLinkingRecurseTime += timeSince(start);
                    //std::cout << "Recursing on self in order to force splitting of endpoints.\n";
                    doLinkingRecurse2(allDetections, allTracklets, searchConfig,
                                      firstEndpoint, secondEndpoint,
                                      newSupportNodes, results, iterationsTillSplit, rangeCache);
                }
                else if (firstEndpointWidth >= secondEndpointWidth) {

                    //"widest" node is first endpoint, recurse twice using its children
                    // in its place.  
                    
                    if ((! firstEndpoint.myTree->hasLeftChild()) && (!firstEndpoint.myTree->hasRightChild())) {
                        throw LSST_EXCEPT(ProgrammerErrorException, "Recursing in a leaf node (first endpoint), must be a bug!");
                    }

                    if (firstEndpoint.myTree->hasLeftChild())
                    {
                        TreeNodeAndTime newTAT(firstEndpoint.myTree->getLeftChild(), firstEndpoint.myTime);
                        //doLinkingRecurseTime += timeSince(start);
                        //std::cout << "Recursing on left child of first endpoint.\n";
                        doLinkingRecurse2(allDetections, allTracklets, searchConfig,
                                          newTAT,secondEndpoint,
                                          newSupportNodes,
                                          results, iterationsTillSplit, rangeCache); 
                        //std::cout << "Returned from recursion on left child of first endpoint.\n";
                    }
                    
                    if (firstEndpoint.myTree->hasRightChild())
                    {
                        TreeNodeAndTime newTAT(firstEndpoint.myTree->getRightChild(), firstEndpoint.myTime);
                        doLinkingRecurseTime += timeSince(start);
                        //std::cout << "recursing on right child of first endpoint..\n";
                        doLinkingRecurse2(allDetections, allTracklets, searchConfig,
                                          newTAT,secondEndpoint,
                                          newSupportNodes,
                                          results, iterationsTillSplit, rangeCache);  
                        //std::cout << "Returned from recursion on right child of first endpoint.\n";
                    }
                }
                else {
                    //"widest" node is second endpoint, recurse twice using its children
                    // in its place
                    
                    if ((!secondEndpoint.myTree->hasLeftChild()) && (!secondEndpoint.myTree->hasRightChild())) {
                        throw LSST_EXCEPT(ProgrammerErrorException, "Recursing in a leaf node (second endpoint), must be a bug!");
                    }

                    if (secondEndpoint.myTree->hasLeftChild())
                    {
                        TreeNodeAndTime newTAT(secondEndpoint.myTree->getLeftChild(), secondEndpoint.myTime);
                        doLinkingRecurseTime += timeSince(start);
                        //std::cout << "Recursing on left child of second endpoint.\n";
                        doLinkingRecurse2(allDetections, allTracklets, searchConfig,
                                          firstEndpoint,newTAT,
                                          newSupportNodes,
                                          results, iterationsTillSplit, rangeCache);
                        //std::cout << "Returned from recursion on left child of second endpoint.\n";
                    }
                    
                    if (secondEndpoint.myTree->hasRightChild())
                    {
                        TreeNodeAndTime newTAT(secondEndpoint.myTree->getRightChild(), secondEndpoint.myTime);
                        doLinkingRecurseTime += timeSince(start);
                        //std::cout << "Recursing on right child of second endpoint.\n";
                        doLinkingRecurse2(allDetections, allTracklets, searchConfig,
                                          firstEndpoint,newTAT,
                                          newSupportNodes,
                                          results, iterationsTillSplit, rangeCache);
                        //std::cout << "Returned from recursion on right child of second endpoint.\n";
                        
                    }
                }
            }                        
        }
    }
}







void debugPrintTimingInfo(const TrackSet &results)
{

    std::cout << " so far we've found " << results.size() << " tracks.\n";
    std::cout << "TIMING STATS: \n---------------------\n";
    std::cout << "getAllDetectionsForTracklet: :\t" << getAllDetectionsForTrackletTime << "sec\n"; 
    std::cout << "getFirstDetectionForTracklet:\t" << getFirstDetectionForTrackletTime << "sec\n"; 
    std::cout << "modifyWithAcceleration:\t" << modifyWithAccelerationTime << "sec\n"; 
    std::cout << "positionAndVelocityRangesOverlapAfterAcceleration:\t" << positionAndVelocityRangesOverlapAfterAccelerationTime << "sec\n"; 
    std::cout << "areCompatible:\t" << areCompatibleTime << "sec\n"; 
    std::cout << "getBestFitVelocityAndAcceleration:\t" << getBestFitVelocityAndAccelerationTime << "sec\n";
    std::cout << "getBestFitVelocityAndAccelerationForTracklets:\t" << getBestFitVelocityAndAccelerationForTrackletsTime << "sec\n"; 
    std::cout << "addBestCompatibleTrackletsAndDetectionsToTrack:\t" << addBestCompatibleTrackletsAndDetectionsToTrackTime << "sec\n"; 
    std::cout << "trackMeetsRequirements:\t" << trackMeetsRequirementsTime << "sec\n"; 
    std::cout << "endpointTrackletsAreCompatible:\t" << endpointTrackletsAreCompatibleTime << "sec\n"; 
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
    getAllDetectionsForTrackletTime = 0; 
    getFirstDetectionForTrackletTime = 0; 
    modifyWithAccelerationTime = 0; 
    positionAndVelocityRangesOverlapAfterAccelerationTime = 0; 
    areCompatibleTime = 0; 
    getBestFitVelocityAndAccelerationTime = 0;
    getBestFitVelocityAndAccelerationForTrackletsTime = 0; 
    addBestCompatibleTrackletsAndDetectionsToTrackTime = 0; 
    trackMeetsRequirementsTime = 0; 
    endpointTrackletsAreCompatibleTime = 0; 
    buildTracksAddToResultsTime = 0; 
    areAllLeavesTime = 0; 
    nodeWidthTime = 0; 
    addAllDetectedObjectsToSetTime = 0; 
    doLinkingRecurseTime = 0;    
}





void doLinking(const std::vector<MopsDetection> &allDetections,
               std::vector<Tracklet> &allTracklets,
               linkTrackletsConfig searchConfig,
               std::map<ImageTime, KDTree <unsigned int> > &trackletTimeToTreeMap,
               TrackSet &results)
{
    /* for every pair of trees, using the set of every intermediate (temporally) tree as a
     * set possible support nodes, call the recursive linker.
     */ 
    /*
     * implementation detail: note that there is a difference between a KDTree
     * and a KDTreeNode.  KDTrees are the "outer layer" and don't really hold
     * data; KDTreeNodes are where the data lives.
     * 
     * In other programs, KDTreeNodes are hidden from the user, but
     * linkTracklets will use them directly.
     *
     * In doLinking_recurse, we will only deal with KDTreeNodes. So in
     * doLinking, we start the searching with the head (i.e. root) node of each
     * tree.
     */
    bool DEBUG = false;

    bool limitedRun = false;
    double limitedRunFirstEndpoint =  49616.273436000003 ;
    double limitedRunSecondEndpoint = 49623.023787999999 ;

    initDebugTimingInfo();

    double iterationTime = std::clock();
    

    if (DEBUG) {
        std:: cout << "all MJDs: ";
        std::map<ImageTime, KDTree<unsigned int> >::const_iterator mapIter;
        for (mapIter = trackletTimeToTreeMap.begin();
             mapIter != trackletTimeToTreeMap.end();
             mapIter++) {
            std::cout << std::setprecision (10) << mapIter->first.getMJD() << " ";
        }
        std::cout << std::endl;
    }

    std::map<ImageTime, KDTree<unsigned int> >::const_iterator firstEndpointIter;
    for (firstEndpointIter = trackletTimeToTreeMap.begin(); 
         firstEndpointIter != trackletTimeToTreeMap.end(); 
         firstEndpointIter++)
    {
        std::map<ImageTime, KDTree<unsigned int> >::const_iterator secondEndpointIter;
        std::map<ImageTime, KDTree<unsigned int> >::const_iterator afterFirstIter = firstEndpointIter;
        afterFirstIter++;

        
        if ((!limitedRun) || (areEqual(firstEndpointIter->first.getMJD(), limitedRunFirstEndpoint))) {
            
            
            for (secondEndpointIter = afterFirstIter; 
                 secondEndpointIter != trackletTimeToTreeMap.end(); 
                 secondEndpointIter++)
            {
                /* if there is sufficient time between the first and second nodes, then try
                   doing the recursive linking using intermediate times as support nodes.
                */
                
                if ((!limitedRun) || 
                    (areEqual(secondEndpointIter->first.getMJD(), limitedRunSecondEndpoint))) {
                    
                    
                    if (secondEndpointIter->first.getMJD() - firstEndpointIter->first.getMJD() 
                        >= searchConfig.minEndpointTimeSeparation) {                
                        
                        if (DEBUG) {
                            time_t rawtime;
                            
                            if (timeSince(iterationTime) > 5) {
                                debugPrintTimingInfo(results);
                            }
                            iterationTime = std::clock();
                            struct tm * timeinfo;
                            
                            time ( &rawtime );
                            timeinfo = localtime ( &rawtime );                    
                            
                            std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;
                            std::cout << "attempting linking between times " << std::setprecision(12)  
                                      << firstEndpointIter->first.getMJD() << " and " 
                                      << std::setprecision(12) << secondEndpointIter->first.getMJD() 
                                      << std::endl;
                        }

                        // get all intermediate points as support nodes.
                
                        /* note that std::maps are sorted by their key, which in this case
                         * is time.  ergo between firstEndpointIter and secondEndpointIter
                         * is *EVERY* tree (and ergo every tracklet) which happened between
                         * the first endpoint's tracklets and the second endpoint's
                         * tracklets.
                         */
                
                        std::vector<TreeNodeAndTime > supportPoints;
                        std::map<ImageTime, KDTree<unsigned int> >::const_iterator supportPointIter;
                        if (DEBUG) {
                            //std::cout << "intermediate times: " ;
                        }
                        for (supportPointIter = afterFirstIter;
                             supportPointIter != secondEndpointIter;
                             supportPointIter++) {
                    
                            /* don't pass along second tracklets which are 'too close' to the endpoints; see
                               linkTracklets.h for more comments */
                            if ((fabs(supportPointIter->first.getMJD() - firstEndpointIter->first.getMJD()) > 
                                 searchConfig.minSupportToEndpointTimeSeparation) 
                                && 
                                (fabs(supportPointIter->first.getMJD() - secondEndpointIter->first.getMJD()) > 
                                 searchConfig.minSupportToEndpointTimeSeparation)) {
                        
                                TreeNodeAndTime tmpTAT(supportPointIter->second.getRootNode(),
                                                       supportPointIter->first);
                                supportPoints.push_back(tmpTAT);
                                if (DEBUG) {
                                    //std::cout << std::setprecision(10) << tmpTAT.myTime << " ";
                                }
                            }
                        }
                
                        if (DEBUG) {std::cout << "\n";}

                        TreeNodeAndTime firstEndpoint(firstEndpointIter->second.getRootNode(), 
                                                      firstEndpointIter->first);
                        TreeNodeAndTime secondEndpoint(secondEndpointIter->second.getRootNode(),
                                                       secondEndpointIter->first);
                
                        //call the recursive linker with the endpoint nodes and support point nodes.

                        // use a cache to avoid doing math in isCompatible!
                        // we use a new cache for each pair of endpoints, because
                        // isCompatible projects the location of the endpoint regions
                        // forward/backwards in time, so there will be 0 reuse between
                        // pairs of endpoint trees.
                        LTCache rangeCache(CACHE_SIZE);
                        //LTCache rangeCache(0);

                        doLinkingRecurse2(allDetections, allTracklets, searchConfig,
                                          firstEndpoint, secondEndpoint,
                                          supportPoints,  
                                          results, ITERATIONS_PER_SPLIT, rangeCache);


                    }
                }
            }
        }
    }
    if (DEBUG) {
        debugPrintTimingInfo(results);
    }
}








TrackSet linkTracklets(const std::vector<MopsDetection> &allDetections,
                       std::vector<Tracklet> &queryTracklets,
                       linkTrackletsConfig searchConfig) {
    TrackSet toRet;
    /*create a sorted list of KDtrees, each tree holding tracklets
      with unique start times (times of first detection in the tracklet).
      
      the points in the trees are in [RA, Dec, RAVelocity, DecVelocity] and the
      returned keys are indices into queryTracklets.
    */
    //std::cout << "setting velocities of tracklets.\n";
    setTrackletVelocities(allDetections, queryTracklets);
    //std::cout << "Building trees on tracklets.\n";
    std::map<ImageTime, KDTree <unsigned int> > trackletTimeToTreeMap;    
    makeTrackletTimeToTreeMap(allDetections, queryTracklets, trackletTimeToTreeMap, searchConfig);
    //std::cout << "Beginning the linking process.\n";
    doLinking(allDetections, queryTracklets, searchConfig, trackletTimeToTreeMap, toRet);
    
    return toRet;
}





}} //close lsst::mops

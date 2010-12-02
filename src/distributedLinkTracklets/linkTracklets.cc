// -*- LSST-C++ -*-
/* jonathan myers */

#include "mpi.h"
#include <ctime>
#include <iomanip>
#include <map>
#include <gsl/gsl_multifit.h>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <sstream>

#include "lsst/mops/daymops/distributedLinkTracklets/linkTracklets.h"
#include "lsst/mops/rmsLineFit.h"
#include "lsst/mops/Exceptions.h"
#include "lsst/mops/KDTree.h"

#include "lsst/mops/daymops/linkTracklets/lruCache.h"


/* since it is expensive to compute the compatible ranges of
 * position/velocity-space associated with each tree node, use a cache which
 * takes a tree node ID and an observation time and reports stores the
 * compatible region.  This seems to give a big performance boost, but you need
 * to set the cache size to something intelligent for your machine/architecture.
 */
#define NODE_COMPATIBILITY_CACHE_SIZE 100000
#define USE_NODE_COMPATIBILITY_CACHE true

/* taking a cue from Kubica, it's only once per ITERATIONS_PER_SPLIT calls to
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

#define POINT_RA           0
#define POINT_DEC          1
#define POINT_RA_VELOCITY   2
#define POINT_DEC_VELOCITY  3

#define uint unsigned int 

namespace lsst {
namespace mops {

/**************************************************
 **************************************************
 **  DISTRIBUTED VARIABLE DEFINITIONS
 **************************************************
 **************************************************/

/**************************************************
 * This struct is used for keeping track of the
 * files written and meta data about them, which
 * will be used when the work is distributed.
 **************************************************/
typedef struct ta{
  int numProcs;
  int sentinel;
}threadArgs;

typedef struct s{
  std::vector<std::pair<uint, uint> > endpoints;
  std::string fileName;
  int numCompatibleEndpoints;
  int numNodes;
  unsigned long long timeUnits;
}workItemFile;

typedef struct treeIdNodeId{
  uint treeId;
  uint nodeId;
}treeIdNodeIdPair;

uint numCachedWorkItems = 0;
unsigned long int totalWorkItems = 0;
unsigned long int totalTimeUnits = 0;
uint currentNumDistributions = 0;
uint numDistributions = 0;
int calculatedTotal = 0;
int globalCacheSize = 0;
pthread_t thread;

#define MAX_WORK_ITEMS 128        /* number of work items per distribution batch */
#define MAX_DISTRIBUTIONS 100000  /* number of work item distribution batches to perform
				     before cleaning up the on-disk work item directory */
#define MPI_COLLECT_TAG 6969      /* MPI tag for track collection */
#define MPI_ANNEAL_TAG 4242       /* MPI tag for annealing trigger */
#define CACHE_SIZE 512
#define WORK_ITEM_CACHE_SIZE 512
#define WRITE_RESULTS_TO_DISK 1   /* workers write tracks to disk or return them */
#define DO_OPRA_1 0               /* perform the optimal page replacement algorithm
				     work item ordering loop or no */

int totalNumTrees, totalNumNodes;
int globalNextWorker = 0;
int numThreads = 0;
std::string rootDir = "workItems";
bool stopAnnealing = false;
bool firstAssignment = true;

/**
 ** END DISTRIBUTED VARS
 **************************************************
 **************************************************/


int doLinkingRecurseVisits = 0, buildTracksAddToResultsVisits = 0, compatibleEndpointsFound = 0;
#define RACE_TO_MAX_COMPATIBLE false
#define MAX_COMPATIBLE_TO_FIND 250000

int rejectedOnVelocity, rejectedOnPosition, wereCompatible;
int cacheHits, cacheMisses;


// THESE ARE FOR DEBUGGING ONLY
#include <ctime>
double getTimeElapsed(clock_t priorEvent)
{
     return ( std::clock() - priorEvent ) / (double)CLOCKS_PER_SEC;
}
double timeSince(clock_t priorEvent)
{
     return ( std::clock() - priorEvent ) / (double)CLOCKS_PER_SEC;
}




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




// a quick macro to remove a lot of ugly typing
#define LTCache LRUCache <GlobalNodeIdAndProjectionId, std::vector <std::vector <double> > >





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
    // NB: since IDs are assigned in order of image time, would it be faster if
    // we did < based on imgId (int rather than double?)
    bool operator< (const ImageTime &other) const {
        return MJD < other.getMJD();
    }

private:
    double MJD;
    uint imgId;
};


class TreeNodeAndTime {
public:
  TreeNodeAndTime(){
  }
    TreeNodeAndTime(KDTreeNode<uint> * tree, ImageTime i) {
    myTree = tree;
    myTime = i;
  }
    KDTreeNode <uint> * myTree;
    ImageTime myTime;
};


/***************************************************************************
 * This is the data that is stored in the node cache.  A node cache entry
 * contains the unique identifiers for the TreeNodeAndTime object (treeId,
 * nodeId) and the TreeNodeAndTime object itself.
 ***************************************************************************/
typedef struct cn{
  
  int treeId, nodeId;
  //TreeNodeAndTime *node;
  TreeNodeAndTime node;

} cachedNode;


/***************************************************
 * From http://www.parashift.com/c++-faq-lite/misc-technical-issues.html#faq-39.1
 ***************************************************/
std::string stringify(double x)
{
  std::ostringstream o;
  o << x;
  return o.str();
} 

/**************************************************
 * From http://www.parashift.com/c++-faq-lite/misc-technical-issues.html#faq-39.2
 **************************************************/
double doublify(std::string s){
  std::istringstream i(s);
  double x;
  i >> x;
  return x;
}


/**************************************************
 * convert string to unsigned int
 **************************************************/
uint intify(std::string s){
  
  std::istringstream i(s);
  uint x;
  i >> x;
  return x;
}


/**************************************************
 * convert string to unsigned long long
 **************************************************/
unsigned long long longlongify(std::string s)
{
  std::istringstream i(s);
  unsigned long long x;
  i >> x;
  return x;
}



void distributeCurrentWorkload(std::vector<std::vector<int> > &assignment,
			       const std::vector<MopsDetection> &allDetections,
			       const std::vector<Tracklet> &allTracklets,
			       linkTrackletsConfig searchConfig,
			       std::map<ImageTime, KDTree <uint> > &trackletTimeToTreeMap,
			       int currentTotalWorkItems, int timeUnits, std::string dir);
  
void *killProcs(threadArgs);

int findFarthestNodeIndex(cachedNode *nodeCache, int cacheSize,
			  const std::vector<std::vector<int> > &finalEndpointOrder);

KDTree<uint> *
findTreeById(const uint id, 
	     std::map<ImageTime, KDTree<uint> > &treeMap);

ImageTime findImageTimeForTree(const uint id, 
			       const std::map<ImageTime, KDTree<uint> > &treeMap);

KDTreeNode<uint> *
getNodeByIDAndTime(KDTree<uint> *myTree, uint nodeId);

int loadNodeFromCache(int treeId, int nodeId, int &cacheSize, cachedNode *nodeCache,
		      std::map<ImageTime, KDTree <uint> > &trackletTimeToTreeMap,
		      const std::vector<std::vector<int> > &finalEndpointOrder,
		      unsigned long long &pageFaults);

/**************************************************
 *
 **************************************************/
std::string timestamp()
{

  time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  std::string retVal(asctime(timeinfo));

  try{
    retVal.erase(retVal.length()-1);
    retVal = "[" + retVal + "]: ";
  }
  catch(std::exception& e){
    std::cerr << "Exception occurred in timestamp().  Output: " << e.what() << std::endl;
  }

  return retVal;
}


/*
 * END
 **************************************************/


// given a tracklet, return the *temporally* earliest detection of this tracklet
// (the one with minimum MJD). 
MopsDetection getFirstDetectionForTracklet(const std::vector<MopsDetection> &allDetections,
                                           const Tracklet &t) 
{
  //double start = std::clock();
    if (t.indices.size() < 1) {
        LSST_EXCEPT(BadParameterException,
                    "linkTracklets::getFirstDetectionForTracklet called with empty tracklet.");
    }

    MopsDetection toRet;
    bool foundOne = false;
    std::set<uint>::const_iterator indexIter;
    for (indexIter = t.indices.begin(); indexIter != t.indices.end(); indexIter++) {
        MopsDetection curDet = allDetections.at(*indexIter);
        if ((!foundOne) || 
            (toRet.getEpochMJD() > curDet.getEpochMJD())) {
            toRet = curDet;
            foundOne = true;
        }
    }

    //getFirstDetectionForTrackletTime += getTimeElapsed(start);
    return toRet;

}





void makeTrackletTimeToTreeMap(const std::vector<MopsDetection> &allDetections,
                               const std::vector<Tracklet> &queryTracklets,
                               std::map<ImageTime, KDTree <uint> > &newMap,
			       int LEAF_SIZE)
{
    newMap.clear();
    //sort all tracklets by their first image time; make PointAndValues from
    //these so we can build a tree.
    // allTrackletPAVsMap will map from image time -> [ all tracklets starting at that image time. ]
    std::map<double, std::vector<PointAndValue<uint> > > allTrackletPAVsMap;
    allTrackletPAVsMap.clear();
    for (uint i = 0; i < queryTracklets.size(); i++) {
        MopsDetection firstDetection =  getFirstDetectionForTracklet(allDetections, queryTracklets.at(i));
        double firstDetectionTime = firstDetection.getEpochMJD();
        PointAndValue<uint> trackletPAV;
        std::vector<double> trackletPoint;

        trackletPoint.push_back(firstDetection.getRA());
        trackletPoint.push_back(firstDetection.getDec());
        trackletPoint.push_back(queryTracklets.at(i).velocityRA);
        trackletPoint.push_back(queryTracklets.at(i).velocityDec);
        trackletPAV.setPoint(trackletPoint);        
	
        trackletPAV.setValue(i);

        allTrackletPAVsMap[firstDetectionTime].push_back(trackletPAV);
	
	//update the number of nodes
	++totalNumNodes;
    }

    // iterate over each time/pointAndValueVec pair and build a corresponding
    // time/KDTree pair.  note that we iterate over a Map which uses MJD as key;
    // Maps sort their data by their key, so we are iterating over all image
    // times in order.
    uint curImageID = 0;

    std::map<double, std::vector<PointAndValue<uint> > >::iterator PAVIter;
    for (PAVIter = allTrackletPAVsMap.begin(); PAVIter != allTrackletPAVsMap.end(); PAVIter++) {
        KDTree<uint> curTree(PAVIter->second, 4, LEAF_SIZE);
        newMap[ImageTime(PAVIter->first, curImageID)] = curTree;
        curImageID++;
	++totalNumTrees;
    }

}




/***************************************************************************
 * Check the cache for the requested node, if it's there, return it, if not
 * load it into the cache using the optimal page replacement algorithm.
 ***************************************************************************/
int loadNodeFromCache(int treeId, int nodeId, int &cacheSize, cachedNode *nodeCache,
		      std::map<ImageTime, KDTree <uint> > &trackletTimeToTreeMap,
		      const std::vector<std::vector<int> > &finalEndpointOrder,
		      unsigned long long &pageFaults)
{
  bool fullCache = true;
  //bool nodeInCache = false;

  //if cache has free slots, load the node in automatically
  if(globalCacheSize < CACHE_SIZE){
    fullCache = false;
  }

  //if this TAT is already in the cache, load it and skip the OPR section
  for(int i=0; i < globalCacheSize; ++i){
    if(nodeCache[i].nodeId == nodeId && nodeCache[i].treeId == treeId){
      //nodeInCache = true;
      return i;
    }
  }

  //std::cerr << timestamp() << "PAGE FAULT" << std::endl;
  ++pageFaults;

  //The TAT was not found in the cach, create new TAT.  Place this new TAT in the
  //cache, evicting the TAT that will be used the farthest in the future
  //if( !nodeInCache ){

  int mostDistantIndex;
  
  //if the cache is full, use OPR to find the index to load this node into
  if(fullCache){
    
    //this is the index in the cache of the node that will be used 
    //farthest in the future
    mostDistantIndex = findFarthestNodeIndex(nodeCache, cacheSize, finalEndpointOrder);
  }
  //if cache is not full, load the node into the next available cache slot
  else{
    mostDistantIndex = globalCacheSize;
    ++globalCacheSize;
  }
  
  //Get the tree and ImageTime for this image (tree) ID
  KDTree<uint> *myTree = findTreeById(treeId, trackletTimeToTreeMap);
  ImageTime it = findImageTimeForTree(treeId, trackletTimeToTreeMap);
  
  //find this node in its tree
  KDTreeNode<uint> *treeNode = NULL;
  if( myTree != NULL ){
    treeNode = getNodeByIDAndTime(myTree, nodeId);
  }
  
  if( treeNode != NULL ){
    nodeCache[mostDistantIndex].nodeId = nodeId;
    nodeCache[mostDistantIndex].treeId = treeId;

    //delete the existing node, if any
    //if(nodeCache[mostDistantIndex].nodeId == -1){
    //delete nodeCache[mostDistantIndex].node;
    //}
    //nodeCache[mostDistantIndex].node = new TreeNodeAndTime(treeNode, it);
    nodeCache[mostDistantIndex].node.myTree = treeNode;
    nodeCache[mostDistantIndex].node.myTime = it;
    
    return mostDistantIndex;
  }
  else{
    std::cerr << timestamp() << "TAT not found, node cache unchanged." << std::endl;
  }

  return (-1);
}






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






std::set<uint> allDetsInTreeNode(KDTreeNode<uint> &t,
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
        // t.getMyData() holds points and values. the "value" part is index into allTracklets.
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








void debugPrint(const TreeNodeAndTime &firstEndpoint, const TreeNodeAndTime &secondEndpoint, 
                std::vector<TreeNodeAndTime> &supportNodes, 
                const std::vector<MopsDetection> &allDetections,
                const std::vector<Tracklet> &allTracklets) 
{
    std::set<uint> leftEndpointDetIds = allDetsInTreeNode(*(firstEndpoint.myTree), 
                                                                  allDetections, allTracklets);
    std::set<uint> rightEndpointDetIds = allDetsInTreeNode(*(secondEndpoint.myTree),
                                                                   allDetections, allTracklets);

    std::cout << "in doLinkingRecurse2,       first endpoint contains     " ;
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




void showNumVisits(KDTreeNode<uint> *tree, uint &totalNodes, uint &totalVisits) 
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
  //double start = std::clock();
    
    detectionsForTracklet.clear();
    std::set<uint>::const_iterator trackletIndexIter;

    for (trackletIndexIter = t.indices.begin();
         trackletIndexIter != t.indices.end();
         trackletIndexIter++) {
        detectionsForTracklet.push_back(allDetections.at(*trackletIndexIter));
    }

    //getAllDetectionsForTrackletTime += getTimeElapsed(start);
}


// modifyWithAcceleration was here
// positionAndVelocityRangesOverlap was here
// extendRangeBackward was here
// areCompatible was here
// setTrackletVelocities was here
// getBestFitVelocityAndAcceleration was here
// getBestFitVelocityAndAccelerationForTracklets was here
// endpointTrackletsAreCompatible was here

class CandidateDetection {
public:
    CandidateDetection() {
        distance = 0; 
        detId = 0;
        parentTrackletId = 0;
    }
    CandidateDetection(double nDistance, uint nDetId, uint nParentTrackletId) {
        distance = nDistance;
        detId = nDetId;
        parentTrackletId = nParentTrackletId;
    }
    double distance;
    uint detId;
    uint parentTrackletId;
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
                                                    const std::vector<uint> candidateTrackletIDs, 
                                                    double RAVelocity, double RAAcceleration, double RAPosition0,
                                                    double DecVelocity, double DecAcceleration, double DecPosition0,
                                                    double time0,
                                                    const linkTrackletsConfig &searchConfig,
                                                    Track &newTrack) 
{
  //double start = std::clock();

    /* the scoreToIDsMap will hold detection MJDs as the key, map from detection
     * time to best-fitting candidate detection at that time
     */
    std::map<double, CandidateDetection > timeToCandidateMap;

    // find the best compatible detection at each unique image time
    std::vector<uint>::const_iterator trackletIDIter;
    std::set<uint>::const_iterator detectionIDIter;
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
            double predRA  = RAPosition0 + RAVelocity*timeOffset + RAAcceleration*timeOffset*timeOffset;
            double predDec = DecPosition0 + DecVelocity*timeOffset + DecAcceleration*timeOffset*timeOffset;
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
    
    /* initialize a list of image times present in the track already. */
    std::set<double> trackMJDs;
    std::set<uint>::const_iterator trackDetectionIndices;
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
    //    addBestCompatibleTrackletsAndDetectionsToTrackTime += timeSince(start);
}







bool trackMeetsRequirements(const std::vector<MopsDetection> & allDetections, 
                            const Track &newTrack, 
                            double RAVelocity, double RAAcceleration, double RAPosition0,
                            double DecVelocity, double DecAcceleration, double DecPosition0,
                            double time0,
                            linkTrackletsConfig searchConfig)
{
  //double start = std::clock();
    bool meetsReqs = true;

    if (newTrack.componentTrackletIndices.size() < searchConfig.minSupportTracklets + 2) {
        meetsReqs = false;
    }

    if (newTrack.componentDetectionIndices.size() < searchConfig.minDetectionsPerTrack) {
        meetsReqs = false;
    }

    //trackMeetsRequirementsTime += timeSince(start);
    return meetsReqs;
    
    
}












/**********************************************************************
 * Matt using so each processor can write out its own results
 **********************************************************************/
void writeMyResults(std::string outFileName, 
		    const std::vector<MopsDetection> &allDets,
		    const std::vector<Tracklet> &allTracklets,
		    TrackSet &toRet, unsigned long long nextNo)
{
  std::ofstream outFile;

  //go through all Tracks
  std::set<Track>::const_iterator curTrack;
  for(curTrack = toRet.componentTracks.begin(); 
      curTrack != toRet.componentTracks.end();
      curTrack++){

    std::string writeName = "results/" + outFileName + "/" + stringify(nextNo) + ".result.txt";
    //std::cerr << timestamp() << "Writing file " << writeName << std::endl;
    outFile.open(writeName.c_str());
    ++nextNo;

    //go through all Track's detections
    std::set<uint>::const_iterator detIter;
    for (detIter = (*curTrack).componentDetectionIndices.begin();
	 detIter != (*curTrack).componentDetectionIndices.end();
	 detIter++) {
      outFile << *detIter << " ";
    }
    outFile << std::endl;
    outFile.close();
  }
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
                             /*TreeNodeAndTime *firstEndpoint,
			       TreeNodeAndTime *secondEndpoint,*/
                             //std::vector<TreeNodeAndTime> &supportNodes,
			     std::vector<treeIdNodeIdPair> &supportNodes,
                             TrackSet & results,
			     int &cacheSize, cachedNode *nodeCache,
			     std::map<ImageTime, KDTree <uint> > &trackletTimeToTreeMap,
			     const std::vector<std::vector<int> > &finalEndpointOrder,
			     unsigned long long &pageFaults)
{
  //double start = std::clock();
  buildTracksAddToResultsVisits++;
  
  if ((firstEndpoint.myTree->isLeaf() == false) ||
      (secondEndpoint.myTree->isLeaf() == false)) {
    LSST_EXCEPT(ProgrammerErrorException, 
		"buildTracksAddToResults got non-leaf nodes, must be a bug!");
  }
  /*for (uint i = 0; i < supportNodes.size(); i++) {
    if (supportNodes.at(i).myTree->isLeaf() == false) {
      LSST_EXCEPT(ProgrammerErrorException, 
		  "buildTracksAddToResults got non-leaf nodes, must be a bug!");            
    }
    }*/
  
  std::vector<PointAndValue<uint> >::const_iterator firstEndpointIter;
  std::vector<PointAndValue<uint> >::const_iterator secondEndpointIter;
  
  /*const std::vector<PointAndValue<uint> > *firstEndpointData = firstEndpoint.myTree->getMyData();
    const std::vector<PointAndValue<uint> > *secondEndpointData = secondEndpoint.myTree->getMyData();*/
  const std::vector<PointAndValue<uint> > *firstEndpointData = firstEndpoint.myTree->getMyData();
  const std::vector<PointAndValue<uint> > *secondEndpointData = secondEndpoint.myTree->getMyData();
  
  for (firstEndpointIter = firstEndpointData->begin();
       firstEndpointIter != firstEndpointData->end();
       firstEndpointIter++) {
    for (secondEndpointIter = secondEndpointData->begin();
	 secondEndpointIter != secondEndpointData->end();
	 secondEndpointIter++) {

      //std::vector<TreeNodeAndTime>::const_iterator supportNodeIter;
      //std::vector<PointAndValue <uint> >::const_iterator supportPointIter;
      std::vector<treeIdNodeIdPair>::const_iterator supportNodeIter;
      std::vector<PointAndValue <uint> >::const_iterator supportPointIter;
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
	
	
	// create a new track with these endpoints
	Track newTrack;
	std::vector<uint> candidateTrackletIDs;
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
	//for (supportNodeIter = supportNodes.begin(); supportNodeIter != supportNodes.end();
	//supportNodeIter++) {
	for (supportNodeIter = supportNodes.begin(); supportNodeIter != supportNodes.end();
	     supportNodeIter++) {
	  treeIdNodeIdPair tini = *supportNodeIter;
	  int i = loadNodeFromCache(tini.treeId, tini.nodeId, cacheSize, 
				    nodeCache, trackletTimeToTreeMap,
				    finalEndpointOrder, pageFaults);
	  
	  const std::vector<PointAndValue <uint> > * curSupportNodeData;
	  //if(!nodeCache[i].node->myTree->isLeaf()){
	  if(!nodeCache[i].node.myTree->isLeaf()){
	    throw LSST_EXCEPT(BadParameterException,
			      std::string(__FUNCTION__) + 
			      std::string(": received non-leaf node as support node."));
	  }
	  //curSupportNodeData = nodeCache[i].node->myTree->getMyData(); 
	  curSupportNodeData = nodeCache[i].node.myTree->getMyData(); 
	  for (supportPointIter  = curSupportNodeData->begin(); 
	       supportPointIter != curSupportNodeData->end();
	       supportPointIter++) {
	    candidateTrackletIDs.push_back(supportPointIter->getValue());
	  }
	}
	
	addBestCompatibleTrackletsAndDetectionsToTrack(allDetections, allTracklets, candidateTrackletIDs,
						       RAVelocity, RAAcceleration, RAPosition0,
						       DecVelocity, DecAcceleration, DecPosition0, time0,
						       searchConfig,
						       newTrack);
	
	if (trackMeetsRequirements(allDetections, newTrack,  
				   RAVelocity, RAAcceleration, RAPosition0, 
				   DecVelocity, DecAcceleration, DecPosition0,
				   time0,searchConfig)) {
	  //std::cerr << timestamp() << "Adding track to results vector..." << std::endl;
	  results.insert(newTrack);
	}
      }
    }    
  }
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







double nodeWidth(KDTreeNode<uint> *node)
{
  double width = 1;
  for (uint i = 0; i < 4; i++) {
    width *= node->getUBounds()->at(i) - node->getLBounds()->at(i);
  }
  return width;    
}




/***************************************************************************
 ***************************************************************************
 **  DISTRIBUTED FUNCTION DEFINITIONS
 ***************************************************************************
 ***************************************************************************/

/***************************************************************************
 * Given an ImageTime id, find the corresponding tree and return it.
 ***************************************************************************/
KDTree<uint> *
findTreeById(const uint id, 
	     std::map<ImageTime, KDTree<uint> > &treeMap)
{
  KDTree<uint> *retVal = NULL;

  std::map<ImageTime, KDTree<uint> >::iterator mapIter;
  for(mapIter = treeMap.begin(); mapIter != treeMap.end(); mapIter++){
    ImageTime it = (*mapIter).first;
    if(it.getImageId() == id){
      retVal = &((*mapIter).second);
      return retVal;
    }
  }

  return retVal;
}


/***************************************************************************
 * Given a ImageTime imageId, find the <ImageTime, KDTree> pair and
 * return the ImageTime object.
 ***************************************************************************/
ImageTime findImageTimeForTree(const uint id, 
			       const std::map<ImageTime, KDTree<uint> > &treeMap)
{

  std::map<ImageTime, KDTree<uint> >::const_iterator mapIter;
  for(mapIter = treeMap.begin(); mapIter != treeMap.end(); mapIter++){
    ImageTime it = (*mapIter).first;
    if(it.getImageId() == id){
      return it;
    }
  }

  ImageTime retVal;
  return retVal;
}



KDTreeNode<uint> *
getNodeByIDAndTime(KDTree<uint> *myTree, uint nodeId)
{
  KDTreeNode<uint> *retVal;

  //search tree for nodeId
  std::vector<KDTreeNode<uint> *> nodeList;
  nodeList.clear();
  
  retVal = myTree->getRootNode();
  nodeList.push_back(retVal);

  while(!nodeList.empty()){

    retVal = nodeList.front();

    //node ID found!
    if(retVal->getId() == nodeId){
      return retVal;
    }

    //add left child
    if(retVal->hasLeftChild()){
      nodeList.push_back(retVal->getLeftChild());
    }

    //add right child
    if(retVal->hasRightChild()){
      nodeList.push_back(retVal->getRightChild());
    }

    //remove element just searched
    nodeList.erase(nodeList.begin());
  }

  return NULL;
}



/**************************************************
 * Remove all work items written to disk
 **************************************************/
void *cleanUp()
{
  std::string cmd = "rm -r workItems/*";
  std::cerr << timestamp() << "Issuing command " << cmd << std::endl;
  int rc = system(cmd.c_str());
  if(rc == 0){
    std::cerr << timestamp() << "disk cleaned" << std::endl;
  }
  else{
    std::cerr << timestamp() << "unable to clean disk" << std::endl;
  }
}



/**************************************************
 * Create a directory for the next round of
 * work items to be written into.
 **************************************************/
void makeNextDirectory()
{
  std::string nextDir = stringify(numDistributions);
  nextDir = rootDir + "/dir_" + nextDir;

  std::string cmd = "mkdir -p " + nextDir;
  std::cerr << timestamp() << "Creating directory " << nextDir << std::endl;
  system(cmd.c_str());
}

/***************************************************************************
 ***************************************************************************
 **  END DISTRIBUTED FUNCTION DEFINITIONS
 ***************************************************************************
 ***************************************************************************/





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
                       //TrackSet & results,
                       int iterationsTillSplit,
                       LTCache &rangeCache, 
		       int numProcs, std::map<ImageTime, KDTree <uint> > &trackletTimeToTreeMap,
		       std::vector<std::vector<int> > &assignment)
{
  if ((RACE_TO_MAX_COMPATIBLE == true) && (compatibleEndpointsFound >= MAX_COMPATIBLE_TO_FIND)) {
    return;
  }
  
  doLinkingRecurseVisits++;
  firstEndpoint.myTree->addVisit();
  
  if (areCompatible(firstEndpoint, secondEndpoint, searchConfig, rangeCache) == false)
    {
      // poor choice of model nodes (endpoint nodes)! give up.
      return;
    }
  else 
    {
      
      std::set<double> uniqueSupportMJDs;
      std::vector<TreeNodeAndTime> newSupportNodes;
      std::vector<TreeNodeAndTime>::iterator supportNodeIter;
      
      /* look through untested support nodes, find the ones that are compatible
       * with the model nodes, add their children to newSupportNodes */
      
      for (supportNodeIter  = supportNodes.begin();
	   supportNodeIter != supportNodes.end();
	   supportNodeIter++) {
	
	if (iterationsTillSplit <= 0) {
	  bool firstEndpointCompatible;
	  bool secondEndpointCompatible;
	  
	  firstEndpointCompatible =  areCompatible(firstEndpoint, *supportNodeIter, searchConfig, rangeCache);
	  secondEndpointCompatible = areCompatible(secondEndpoint, *supportNodeIter, searchConfig, rangeCache);
          
	  if (firstEndpointCompatible && secondEndpointCompatible)
	    {   
	      if (supportNodeIter->myTree->isLeaf()) {
		newSupportNodes.push_back(*supportNodeIter);
		uniqueSupportMJDs.insert(supportNodeIter->myTime.getMJD());
	      }
	      else {
		if (supportNodeIter->myTree->hasLeftChild()) {
		  
		  TreeNodeAndTime newTAT(supportNodeIter->myTree->getLeftChild(), supportNodeIter->myTime);
		  /*
		    uint treeId = supportNodeIter->myTime.getImageId();
		    uint parentId = supportNodeIter->myTree->getId();
		    uint childId = newTAT.myTree->getId();
		    parentTable[treeId][childId] = parentId;
		  */
		  newSupportNodes.push_back(newTAT);
		  uniqueSupportMJDs.insert(supportNodeIter->myTime.getMJD());
		}
		if (supportNodeIter->myTree->hasRightChild()) {
		  
		  TreeNodeAndTime newTAT(supportNodeIter->myTree->getRightChild(), supportNodeIter->myTime);
		  /*
		    uint treeId = firstEndpoint.myTime.getImageId();
		    uint parentId = firstEndpoint.myTree->getId();
		    uint childId = newTAT.myTree->getId();
		    parentTable[treeId][childId] = parentId;
		  */
		  newSupportNodes.push_back(newTAT);
		  uniqueSupportMJDs.insert(supportNodeIter->myTime.getMJD());
		}
	      }
	    }
	}
	else {
	  // iterationsTillSplit is non-zero; 
	  // don't do any real work; just copy the previous support nodes. we'll do this
	  // momentarily; but go ahead and add their support MJDs now.
	  uniqueSupportMJDs.insert(supportNodeIter->myTime.getMJD());                
	}
      }
      
      
      if (iterationsTillSplit <= 0) {
	iterationsTillSplit = ITERATIONS_PER_SPLIT;
      }
      else{
	// we still need to get newSupportNodes set up. just use the old ones.
	newSupportNodes = supportNodes;
      }
      
      if (uniqueSupportMJDs.size() < searchConfig.minSupportTracklets) {
	// we can't possibly have enough distinct support tracklets between endpoints
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
	    
	    /**************************************************************
	     **************************************************************
	     **  DISTRIBUTED VERSION DOES NOT CALL buildTracksAddTResults
	     **  DIRECTLY, IT CREATES A SET OF WORK ITEMS AND THEN 
	     **  DISTRIBUTES THEM AFTER ANALYSIS
	     **************************************************************
	     **************************************************************/
	    
	    /************************************************************
	     * Calculate the number of compatible endpoints in this set.
	     * This allows us to predict the amount of work required by
	     * this work item and distributed it accordingly.
	     ************************************************************/
	    double RAVelocity, DecVelocity, RAAcceleration, DecAcceleration;
	    double RAPosition0, DecPosition0;
	    double time0;
	    int numCompatible = 0;
	    
	    std::vector<PointAndValue<uint> >::const_iterator firstEndpointIter;
	    std::vector<PointAndValue<uint> >::const_iterator secondEndpointIter;
	    
	    for (firstEndpointIter = firstEndpoint.myTree->getMyData()->begin();
		 firstEndpointIter != firstEndpoint.myTree->getMyData()->end();
		 firstEndpointIter++) {
	      
	      for (secondEndpointIter = secondEndpoint.myTree->getMyData()->begin();
		   secondEndpointIter != secondEndpoint.myTree->getMyData()->end();
		   secondEndpointIter++) {
		
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
		  ++numCompatible;
		  ++compatibleEndpointsFound;
		}
	      }
	    }
	    
	    /****************************************************
	     * Write the information of the endpoints and their 
	     * support points to file, this is a work item
	     * as defined in the terminology.
	     ****************************************************/
	    std::vector<std::pair<uint, uint> > _endpoints;
	    _endpoints.clear();
	    
	    //first item is first endpoint
	    uint treeId = firstEndpoint.myTime.getImageId();
	    uint nodeId = firstEndpoint.myTree->getId();
	    _endpoints.push_back(std::make_pair(treeId, nodeId));
	    
	    //second item is second endpoint
	    treeId = secondEndpoint.myTime.getImageId();
	    nodeId = secondEndpoint.myTree->getId();
	    _endpoints.push_back(std::make_pair(treeId, nodeId));
	    
	    //add all endpoints
	    for(uint m=0; m < newSupportNodes.size(); ++m){
	      treeId = newSupportNodes.at(m).myTime.getImageId();
	      nodeId = newSupportNodes.at(m).myTree->getId();
	      _endpoints.push_back(std::make_pair(treeId, nodeId));
	    }
	    
	    //create workItemFile struct to go on list
	    workItemFile myItem;
	    myItem.endpoints = _endpoints;
	    myItem.numNodes = 2 + newSupportNodes.size();
	    myItem.numCompatibleEndpoints = numCompatible;
	    myItem.timeUnits = myItem.numNodes * myItem.numCompatibleEndpoints;
	    totalTimeUnits += myItem.timeUnits;
	    
	    //write work item details to file
	    std::string dirNo = stringify(numDistributions);
	    dirNo = rootDir + "/dir_" + dirNo;

	    std::ofstream outFile;
	    std::string outFileName = stringify((double)totalWorkItems);
	    outFileName = dirNo + "/" + outFileName;
	    outFile.open(outFileName.c_str());
	    
	    for (uint i = 0; i < myItem.endpoints.size(); i++) {
	      outFile << myItem.endpoints.at(i).first << " " << myItem.endpoints.at(i).second << " ";
	    }
	    outFile << "; " << myItem.numNodes << "; " << myItem.numCompatibleEndpoints << "; " << myItem.timeUnits;
	    outFile.close();
	    
	    //update the work item set assignment vectors
	    assignment.at(globalNextWorker).push_back(totalWorkItems);
	    
	    ++globalNextWorker;
	    if(globalNextWorker >= (numProcs-1)){
	      globalNextWorker = 0;
	    }
	    
	    ++totalWorkItems;
	    
	    //distribute workload periodically
	    if(totalWorkItems >= MAX_WORK_ITEMS){
	      
	      distributeCurrentWorkload(assignment, allDetections, allTracklets, searchConfig, 
					trackletTimeToTreeMap, totalWorkItems,
					totalTimeUnits, dirNo);

	      //update and reset assignment-specific variables
	      ++numDistributions;
	      ++currentNumDistributions;
	      totalWorkItems = 0;
	      totalTimeUnits = 0;

	      //create the directory for the next set of work items to be stored
	      makeNextDirectory();
	      
	      //clear our our work item cache and start afresh
	      for(uint i=0; i < assignment.size(); ++i){
		assignment.at(i).clear();
	      }
	      
	      //tell the workers that this round of work item distribution is complete and
	      //collect any tracks that were generated
	      threadArgs ta;
	      ta.numProcs = numProcs;
	      ta.sentinel = -2;
	      killProcs(ta);

	      //clean up work item files, if necessary
	      if(currentNumDistributions > MAX_DISTRIBUTIONS){
		std::cout << timestamp() << "Cleaning the work item file directory..." << std::endl;
		cleanUp();
		currentNumDistributions = 0;
	      }
	    }
	    
	    if ((RACE_TO_MAX_COMPATIBLE == true) && (compatibleEndpointsFound >= MAX_COMPATIBLE_TO_FIND)) {
	      return;
	    }
	    
	    /**************************************************
	     **************************************************
	     ** This is where buildTracksAddToResults is called
	     ** in the sequential program. -- END DISTRIBUTED
	     **************************************************
	     **************************************************/
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
		  doLinkingRecurse2(allDetections, allTracklets, searchConfig,
				    firstEndpoint,secondEndpoint,
				    newSupportNodes,
				    iterationsTillSplit, rangeCache, 
				    numProcs, trackletTimeToTreeMap, assignment);
		  
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
			      
		      doLinkingRecurse2(allDetections, allTracklets, searchConfig,
					newTAT, secondEndpoint,
					newSupportNodes,
					iterationsTillSplit, rangeCache, 
					numProcs, trackletTimeToTreeMap, assignment);
                    }
		  
		  if (firstEndpoint.myTree->hasRightChild())
                    {
		      TreeNodeAndTime newTAT(firstEndpoint.myTree->getRightChild(), firstEndpoint.myTime);
			      
		      doLinkingRecurse2(allDetections, allTracklets, searchConfig,
					newTAT, secondEndpoint,
					newSupportNodes,
					iterationsTillSplit, rangeCache, 
					numProcs, trackletTimeToTreeMap, assignment);
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

		      doLinkingRecurse2(allDetections, allTracklets, searchConfig,
					firstEndpoint,newTAT,
					newSupportNodes,
					iterationsTillSplit, rangeCache, 
					numProcs, trackletTimeToTreeMap, assignment);
                    }
		  
		  if (secondEndpoint.myTree->hasRightChild())
                    {
		      TreeNodeAndTime newTAT(secondEndpoint.myTree->getRightChild(), secondEndpoint.myTime);
		      
		      doLinkingRecurse2(allDetections, allTracklets, searchConfig,
					firstEndpoint,newTAT,
					newSupportNodes,
					iterationsTillSplit, rangeCache, 
					numProcs, trackletTimeToTreeMap, assignment);
                    }
                }
            }                        
        }
    }
}



//doLinking was here






/***************************************************************************
 ***************************************************************************
 **   BEGIN DISTRIBUTED LINKTRACKLETS FUNCTION DEFINITIONS
 **     --MATT CLEVELAND
 ***************************************************************************
 ***************************************************************************/


/**************************************************
 * collect all Tracks generated by slaves
 **************************************************/
//TrackSet collectTrackResults(int numProcs)
void collectTrackResults(long numProcs)
{
  MPI_Status status;

  int _from;
  for(_from = 1; _from < numProcs; ++_from){

    //determine number of tracks to expect
    int numTracks;
    MPI_Recv(&numTracks, 1, MPI_INT, _from, MPI_COLLECT_TAG, MPI_COMM_WORLD, &status);
    std::cerr << timestamp() << "Expecting " << numTracks << " tracks from the worker " << _from << "." << std::endl;
    
    for(int i=0; i < numTracks; ++i){
      
      //receive tracklet indices
      int numTracklets;
      MPI_Recv(&numTracklets, 1, MPI_INT, _from, MPI_COLLECT_TAG, MPI_COMM_WORLD, &status);
      
      int trackletIndices[numTracklets];
      MPI_Recv(&trackletIndices, numTracklets, MPI_INT, _from, MPI_COLLECT_TAG, MPI_COMM_WORLD, &status);
      
      std::set<uint> tSet;
      for(int k=0; k < numTracklets; ++k){
	tSet.insert(trackletIndices[k]);
      }
      
      //receive detection indices
      int numDetections;
      MPI_Recv(&numDetections, 1, MPI_INT, _from, MPI_COLLECT_TAG, MPI_COMM_WORLD, &status);    
      
      int detectionIndices[numDetections];
      MPI_Recv(&detectionIndices, numDetections, MPI_INT, _from, MPI_COLLECT_TAG, MPI_COMM_WORLD, &status);
      
      std::set<uint> dSet;
      for(int k=0; k < numDetections; ++k){
	dSet.insert(detectionIndices[k]);
      }
      
      //add the new track to our return set
      Track *thisTrack = new Track();
      thisTrack->componentTrackletIndices = tSet;
      thisTrack->componentDetectionIndices = dSet;
      
      toRet.insert((*thisTrack));
    }
  }
}



/***************************************************************************
 * @@Helper function for simulated annealing@@
 *
 * 
 ***************************************************************************/
workItemFile extractWorkItemFromString(std::string data)
{
  workItemFile retVal;

  std::string token;
  std::istringstream iss(data);
  std::vector<std::string> allTokens;

  while(getline(iss, token, ';')){
    allTokens.push_back(token);
  }

  //ensure valid string format was supplied
  if(allTokens.size() < 4){
    retVal.numNodes = -1;
    return retVal;
  }
    
  std::string pt;
  std::istringstream iss2(allTokens.at(0));

  //load all endpoints
  int i = 0;
  uint ti;
  uint ni;
  while(getline(iss2, pt, ' ')){
    if(i % 2 == 0){
      ti = intify(pt);
    }
    else{
      ni = intify(pt);
      retVal.endpoints.push_back(std::make_pair(ti, ni));
    }

    ++i;
  }

  //load numNodes
  retVal.numNodes = intify(allTokens.at(1));
  
  //load numCompatibleEndpoints
  retVal.numCompatibleEndpoints = intify(allTokens.at(2));

  //load time units
  retVal.timeUnits = longlongify(allTokens.at(3));

  return retVal; 
}



/***************************************************************************
 * @@Helper function for simulated annealing@@
 *
 * 
 ***************************************************************************/
workItemFile readWorkItemFileFromDisk(int fileNo, workItemFile *workItemCache, 
				      std::string dir)
{
  //open file and read its contents
  workItemFile retVal;

  std::string fileName = stringify((double)fileNo);
  fileName = dir + "/" + fileName;

  //check to see if this work item is in the cache
  if(numCachedWorkItems > 0){
    for(uint i=0; i < numCachedWorkItems; ++i){
      if(fileName == workItemCache[i].fileName){
	retVal.endpoints = workItemCache[i].endpoints;
	retVal.numNodes = workItemCache[i].numNodes;
	retVal.numCompatibleEndpoints = workItemCache[i].numCompatibleEndpoints;
	retVal.timeUnits = workItemCache[i].timeUnits;
	return retVal;
      }
    }
  }

  std::ifstream myFile;
  myFile.open(fileName.c_str());

  std::string line;
  std::vector<std::string> lines;
  lines.clear();

  if (myFile.is_open()){
    while (!myFile.eof()){
      getline (myFile,line);
      lines.push_back(line);
    }
    myFile.close();
  }
  
  //extract work item data from line
  if(lines.size() == 0){
    std::cerr << timestamp() << "Empty work item file encountered." << std::endl;
  }
  else if(lines.size() == 1){
    retVal = extractWorkItemFromString(lines.at(0));

    //update workFileCache
    int nodeIndex = numCachedWorkItems;
    if(numCachedWorkItems >= CACHE_SIZE){
      nodeIndex = rand() % CACHE_SIZE;
    }
    else{
      ++numCachedWorkItems;
    }

    workItemCache[nodeIndex].endpoints = retVal.endpoints;
    workItemCache[nodeIndex].numNodes = retVal.numNodes;
    workItemCache[nodeIndex].numCompatibleEndpoints = retVal.numCompatibleEndpoints;
    workItemCache[nodeIndex].timeUnits = retVal.timeUnits;
    workItemCache[nodeIndex].fileName = fileName;
  }
  else{
    std::cerr << timestamp() << "Corrupted work item file encountered." << std::endl;
  }

  return retVal;
}



/***************************************************************************
 * @@Helper function for simulated annealing@@
 *
 * Determine the cost of this set assignment using the exeuction model.
 ***************************************************************************/
unsigned long long executionModelCost(int itemNo, workItemFile *workItemCache,
				      std::string dir)
{
  //load the work item from it's file on disk <itemNo>
  workItemFile myItem = readWorkItemFileFromDisk(itemNo, workItemCache, dir);

  //initialize to computational model cost
  unsigned long long cost = myItem.timeUnits;
  
  //add in communication cost
  cost += (myItem.endpoints.size() * 2);

  return cost;
}



/***************************************************************************
 * @@Helper function for simulated annealing@@
 *
 * Take the "assignment" and shuffle the work items it contains to create
 * a "neighbor" of this distribution.
 ***************************************************************************/
void generateNeighbor(std::vector<std::vector<int> > &assignment, 
		      std::vector<std::vector<int> > &retVal,
		      int currentTimeUnits,
		      workItemFile *workItemCache, 
		      std::string dir)
{
  std::vector<std::vector<int> > mixedSet;
  //std::cerr << timestamp() << "Entering generateNeighbor..." << std::endl;

  //copy assignment into mixedSet so we can muck with it
  uint numWorkers = assignment.size();
  for(uint i=0; i < numWorkers; ++i){
    mixedSet.push_back(assignment.at(i));
  }

  //shuffle all the work items assigned to each worker
  std::random_shuffle(mixedSet.begin(), mixedSet.end());
  for(uint i=0; i < mixedSet.size(); ++i){
    std::random_shuffle(mixedSet.at(i).begin(), mixedSet.at(i).end());
  }

  unsigned long long timePerWorker = (currentTimeUnits / numWorkers);
  int nextWorkItem = 0;
  uint workItemNodeSize = mixedSet.at(nextWorkItem).size();
  uint numAssignedWorkItems = 0;
  uint workItemNodeIndex = 0;
  unsigned long long assignedTime = 0;

  //each worker gets their share of work items
  for(uint i=0; i < numWorkers; ++i){
    
    std::vector<int> myItems;

    //keep giving this worker work until they've gotten 
    //a balanced amount of the total time
    while(assignedTime < timePerWorker){
      
      if(workItemNodeIndex >= workItemNodeSize){
	++nextWorkItem;

	if(nextWorkItem >= mixedSet.size()){
	  break;
	}

	workItemNodeIndex = 0;
	workItemNodeSize = mixedSet.at(nextWorkItem).size();
      }
      
      if(workItemNodeIndex < mixedSet.at(nextWorkItem).size()){
	myItems.push_back(mixedSet.at(nextWorkItem).at(workItemNodeIndex));
	assignedTime += executionModelCost(mixedSet.at(nextWorkItem).at(workItemNodeIndex), workItemCache, dir);
	++workItemNodeIndex;
	++numAssignedWorkItems;
      }
    }

    assignedTime = 0;
    
    retVal.push_back(myItems);
  }

  //assign any leftover work items
  uint numRemaining = totalWorkItems - numAssignedWorkItems;

  if(numRemaining > 0){

    int nextWorker = 0;
    int n = mixedSet.size() - 1;
    int remainder = numRemaining - mixedSet.at(n).size();

    while(remainder > 0){
      --n;
      remainder = remainder - mixedSet.at(n).size();
    }
    
    //this is the index into the nth work item set we should start with
    if(n != (mixedSet.size() - 1)){
      remainder = -remainder;
    }
    else{
      remainder = mixedSet.at(n).size() - numRemaining;
    }

    for(uint i = n; i < mixedSet.size(); ++i){
      for(uint j = remainder; j < mixedSet.at(i).size(); ++j){
	retVal.at(nextWorker).push_back(mixedSet.at(i).at(j));
	
	//assign the remaining work items cyclically
	++nextWorker;
	++numAssignedWorkItems;
	if(nextWorker == numWorkers){
	  nextWorker = 0;
	}
      }
      remainder = 0;
    }
  }
}



/***************************************************************************
 * @@Helper function for simulated annealing@@
 *
 * Determine the cost of this set assignment using the exeuction model.
 ***************************************************************************/
unsigned long long calculateSetAssignmentCost(std::vector<std::vector<int> > &assignment,
					      workItemFile *workItemCache, 
					      std::string dir)
{

  //L is the constant value representing the cost to load a node from memory
  uint L = 1000;

  //this is the max cost from this work item set assignment
  unsigned long long maxCost = 0;

  for(uint i=0; i < assignment.size(); ++i){
    
    unsigned long long cost = 0;

    //look through all the nodes of the adjacent work item and find number of shared nodes
    for(uint j=0; j < assignment.at(i).size(); ++j){
      cost += executionModelCost(assignment.at(i).at(j), workItemCache, dir);
      
      workItemFile oneTemp = readWorkItemFileFromDisk(assignment.at(i).at(j), workItemCache, dir);
      uint matches = 0;
      uint numNodes = oneTemp.numNodes;;
      
      // memory (IO) model cost (this only goes here because the IO
      // cost varies at the WORK ITEM SET ASSIGNMENT level, not the work
      // item level)
      if(j == 0){
	cost+=(numNodes * L);
      }
      else{

	//think: should this be (j+1)?
	workItemFile twoTemp = readWorkItemFileFromDisk(assignment.at(i).at((j-1)), workItemCache, dir);

	for(uint k=0; k <  oneTemp.endpoints.size(); ++k){
	  for(uint l=0; l < twoTemp.endpoints.size(); ++l){
	    if((oneTemp.endpoints.at(k).first == twoTemp.endpoints.at(l).first) &&
	       (oneTemp.endpoints.at(k).second == twoTemp.endpoints.at(l).second)){
	      ++matches;
	      break;
	    }
	  }
	}
     
	cost+=(L * (numNodes - matches));

      }/*end else*/
    }/* end j */

    //keep track of the greatest cost encountered
    if(cost > maxCost){
      maxCost = cost;
    }
    
  }/* end i */

  return maxCost;
}


/***************************************************************************
 * @@Helper function for simulated annealing@@
 *
 * Generate a random number between 0 and 1
 * return a uniform number in [0,1].
 *
 * source: http://www.cplusplus.com/forum/beginner/7445/
 ***************************************************************************/
double unifRand()
{
  return rand() / double(RAND_MAX);
}


/***************************************************************************
 * @@Helper function for simulated annealing@@
 *
 * Generate a random number in a real interval.
 * param a one end point of the interval
 * param b the other end of the interval
 * return a inform rand numberin [a,b].
 *
 * source: http://www.cplusplus.com/forum/beginner/7445/
 ***************************************************************************/
double unifRand(double a, double b)
{
  return (b-a)*unifRand() + a;
}



/***************************************************************************
 * @@Helper function for simulated annealing@@
 *
 * Consider the cost of "cost" as compared to the "bestCost".  Weight its 
 * using a random weighting.  If the weighted value is greater than bestCost
 * return true, otherwise false.
 ***************************************************************************/
bool updateCurrent(uint cost, uint bestCost, uint T)
{
  double randomWeight;
  switch(T){
  case 10:
    randomWeight = unifRand(.3, 1);
    break;
  case 9:
    randomWeight = unifRand(.325, 1);
    break;
  case 8:
    randomWeight = unifRand(.375, 1);
    break;
  case 7:
    randomWeight = unifRand(.45, 1);
    break;
  case 6:
    randomWeight = unifRand(.55, 1);
    break;
  case 5:
    randomWeight = unifRand(.675, 1);
    break;
  case 4:
    randomWeight = unifRand(.825, 1);
    break;
  case 3:
    randomWeight = unifRand(.9, 1);
    break;
  case 2:
    randomWeight = unifRand(.925, 1);
    break;
  case 1:
    randomWeight = unifRand(.95, 1);
    break;
  default:
    randomWeight = unifRand(.975, 1);
    break;
  }

  double _tempCost = (double) cost;
  double _tempBest = (double) bestCost;

  /*if(((_tempCost * randomWeight) < _tempBest) ||
    (cost < bestCost)){*/
  if((_tempCost * randomWeight) < _tempBest){
    return true;
  }
  else{
    return false;
  }
}



/***************************************************************************
 * This algorithm is an implementation of simulated annealing.  The seed set
 * is the naive distribution of work items, that is every work item is 
 * assigned to the next worker node in a round-robin manner.  From this 
 * seed, the work item set assignment is shuffled and the cost of the 
 * assignment is calculated.  In general, the work item set assignment with
 * the lowest cost is kept as the `current' work item set assignment, however
 * according to a certain random factor, some non-optimal work item set
 * assignments will be selected as `current'.  This is the key to simulated
 * annealing and helps prevent the algorithm from settling on a local
 * optimum instead of a global one.
 ***************************************************************************/
void determineWorkItemSets(std::vector<std::vector<int> > &assignment, 
			   int currentTotalWorkItems, int currentTimeUnits, 
			   workItemFile *workItemCache, std::string dir)
{
  //one of the processors is the master
  int numWorkers = assignment.size();

  //Vector (per worker) of vectors of work items.  This is a work item
  //set assignment
  std::vector<std::vector<int> > currentAss;
  for(uint i=0; i < assignment.size(); ++i){
    currentAss.push_back(assignment.at(i));
  }

  //number of iterations to consider different work item set 
  //assignment permutations
  int maxIterations = 10, throttle = 500;
  int iteration = 0;
  uint T = 10;

  //record-keeping variables for the simulated annealing algorithm
  unsigned long long bestCost = ULONG_MAX;

  //ensure at least one round of "intelligent" work item 
  //distribution is considered
  bool singleRun = false;
  if(stopAnnealing){
    singleRun = true;
  }

  while(!stopAnnealing || singleRun){

    //this is the work item set assignment created as a result of 
    //shuffling `current'
    std::vector<std::vector<int> > trialAss;
    generateNeighbor(currentAss, trialAss, currentTimeUnits, workItemCache, dir);
    unsigned long long trialCost = calculateSetAssignmentCost(trialAss, workItemCache, dir);

    //this is the new best set assignment
    if( trialCost < bestCost ){
      std::cerr << timestamp() << "updating best work item set assignment in SA. Trial cost is " 
		<< trialCost << " bestCost is " << bestCost << std::endl;
      bestCost = trialCost;
      
      //clear the current best 
      for(uint i=0; i < assignment.size(); ++i){
	assignment.at(i).clear();
      }
      assignment.clear();
      
      //replace current best with the new best
      for(uint i=0; i < trialAss.size(); ++i){
	assignment.push_back(trialAss.at(i));
      }
    }

    //check if we should update the current set assignment
    if(updateCurrent(trialCost, bestCost, T)){

      //clear the current set assignment
      for(uint i=0; i < currentAss.size(); ++i){
	currentAss.at(i).clear();
      }
      currentAss.clear();
      
      //replace current set assignment with trial set assignment
      for(uint i=0; i < trialAss.size(); ++i){
	currentAss.push_back(trialAss.at(i));
      }
    }
    
    //this catch is for the first work item set assignment, there will
    //be no thread to trigger the "stopAnnealing" loop variable yet
    if((firstAssignment && iteration > maxIterations)){
      firstAssignment = false;
      break;
    }

    //this means that the next distribution is ready to go out, kick out
    //of the loop
    if(singleRun){
      singleRun = false;
      break;
    }

    //"temperature" for P(e, e', T) function, never < 0
    if(T >= 1)
      --T;

    //iterations
    ++iteration;
  }
  
  //reset the annealing sentinel
  stopAnnealing = false;
}



/***************************************************************************
 * Wait to hear from all workers that they have completed processing their
 * work items.  When have heard from all, tell simulated annealing algorithm
 * to stop.
 ****************************************************************************/
void *annealingSentinel(void *arg)
{
  int numProcs = (int)arg;
  MPI_Status status;

  for(int i=1; i <= numProcs; ++i){
    int annealBuff;
    std::cerr << timestamp() << "pthread MPI loop " << i << "..." << std::endl;
    MPI_Recv(&annealBuff, 1, MPI_INT, i, MPI_ANNEAL_TAG, MPI_COMM_WORLD, &status);
    std::cerr << timestamp() << "pthread MPI loop " << i << " has received the go-ahead" << std::endl;
  }

  stopAnnealing = true;
  pthread_detach(pthread_self());
}



/************************************************************
 * Distribute work items to worker nodes according to the 
 * distribution defined in "workItemSets".
 ************************************************************/
void distributeWorkload(std::vector<std::vector<int> > &workItemSets,
			workItemFile *workItemCache, std::string dir)
{
  uint numProcs = workItemSets.size();

  //std::cerr << timestamp() << "Master distributing " << totalWorkItems << " work items to " << numProcs << " workers."  <<  std::endl;

  long int workUnitsList[numProcs];
  for(int i=0; i<numProcs; ++i){
    workUnitsList[i] = 0;
  }

  for(uint i=0; i < workItemSets.size(); ++i){

    //corner case where a worker wasn't assigned any work items, inform him that he 
    //should expect zero items
    if(workItemSets.at(i).size() == 0){
      int numEndpoints = 0;
      MPI_Send(&numEndpoints, 1, MPI_INT, (i+1), i, MPI_COMM_WORLD);
    }
    else{
      for(uint j=0; j < workItemSets.at(i).size(); ++j){
	
	workItemFile workItem = readWorkItemFileFromDisk(workItemSets.at(i).at(j), workItemCache, dir);
	
	//send number of endpoints to recipient
	int numEndpoints = workItem.endpoints.size();

	//send to (i+1), i starts at 0 but proc 0 doesn't process work items
	MPI_Send(&numEndpoints, 1, MPI_INT, (i+1), j, MPI_COMM_WORLD);
	
	//keep track of predicted time required at each processor
	workUnitsList[i] += workItem.timeUnits;
	
	//aggregate data to send from each endpoint
	int bufferSize = numEndpoints * 2;
	int sendBuffer[bufferSize];
	int index = 0;
	
	for(int k=0; k < numEndpoints; ++k){
	  
	  //collect this endpoint's time index, this allows for proper
	  //tree identification
	  sendBuffer[index] = workItem.endpoints.at(k).first;
	  ++index;
	  
	  //collect this endpoint's ID so we can find this node
	  //in the tree specified by timeIndex
	  sendBuffer[index] = workItem.endpoints.at(k).second;
	  ++index;
	} 
	
	//send to (i+1), i starts at 0 but proc 0 doesn't process work items
	MPI_Send(&sendBuffer, bufferSize, MPI_INT, (i+1), j, MPI_COMM_WORLD);
      }
    }
  } /* end file transmission */


  //just for debug
  //std::cerr << timestamp() << "Master distributed workload as such:" << std::endl;
  for(int i=1; i <= numProcs; ++i){
    std::string myTime("Processor ");
    myTime = timestamp() + myTime + stringify(i) + " got " + stringify(workUnitsList[(i-1)]) 
      + " time units and " + stringify(workItemSets.at((i-1)).size()) + " total work items.\n";
    std::cerr << myTime;
  }

  //this thread waits to hear back from the workers that they have completed
  //processing on the work items they were assigned
  int rc = pthread_create(&thread, NULL, annealingSentinel, (void*)numProcs);
  if(rc){
    std::cerr << "pthread_create return code is" << rc << ", aborting." << std::endl;
    exit(rc);
  }
}



/**********************************************************************
 * Tell all the slave nodes they can stop working.
 **********************************************************************/
//void killProcs(int numProcs, int sentinel, TrackSet &TrackResults){
//void killProcs(int numProcs, int sentinel){
void *killProcs(threadArgs ta)
{
  int _sentinel = ta.sentinel;
  
  for(uint i=1; i < ta.numProcs; ++i){
    std::cerr << timestamp() << "sending kill sentinel to " << i << std::endl;
    MPI_Send(&_sentinel, 1, MPI_INT, i, i, MPI_COMM_WORLD);
  }

  //collect results and add to return vector if workers not writing to disk
  if(!WRITE_RESULTS_TO_DISK){
    collectTrackResults(ta.numProcs);
  }
}



/***************************************************************************
 * Look through the cache and see which cached node is used farthest in
 * the future.  If a cached node is found that is never used again, there
 * is no need to check any more nodes.
 ***************************************************************************/
int findFarthestNodeIndex(cachedNode *nodeCache, int cacheSize,
			  const std::vector<std::vector<int> > &finalEndpointOrder)
{
  int farthest = -1;

  std::vector<std::vector<int> >::const_iterator endpointIter;

  //check the cache to see if this node exists
  for(int i=0; i < globalCacheSize; ++i){
    
    bool found = false;

    //look at each work item in order
    for(uint j=0; j< finalEndpointOrder.size(); ++j){
      
      //go through this work item's nodes
      for(uint k=0; k < finalEndpointOrder.at(j).size() && !found; ++k){

	if(finalEndpointOrder.at(j).at(k) == nodeCache[i].treeId && 
	   finalEndpointOrder.at(j).at((k+1)) == nodeCache[i].nodeId){

	  int nodeOffset = j + k;
	  if(nodeOffset > farthest){
	    farthest = i;
	  }
	  found = true;
	}

	//treeId/nodeId pair increment
	++k;
      }
    }
    
    //this means the node is never used again, evict it and move on
    if(!found){
      farthest = i;
      break;
    }
  }

  return farthest;
}






/************************************************************
 * Slave processors wait for signals from master.
 * They accept the data and process it in doLinkingRecurse2.
 ************************************************************/
void waitForTask(int rank,
		 const std::vector<MopsDetection> &allDetections,
		 std::vector<Tracklet> &allTracklets, 
		 linkTrackletsConfig searchConfig)
{
  //this is the set of tracks this function returns
  TrackSet myTracks;

  //this is THE node cache, declare and initialize it
  cachedNode nodeCache[CACHE_SIZE];
  for(uint i=0; i < CACHE_SIZE; ++i){
    nodeCache[i].nodeId = -1;
    nodeCache[i].treeId = -1;
  }

  //number of elements in the cache
  int cacheSize = 0;

  //Receive LNS size from master, this tells us with which leaf node size
  //we should build our tree
  int LEAF_SIZE = 32;
  /*  
      MPI_Status leafStatus;
      MPI_Recv(&LEAF_SIZE, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &leafStatus);
      std::cerr << "Worker received LNS assignment " << LEAF_SIZE << std::endl;
  */
  
  //need to make the time/tree map so we can identify common Node IDs and 
  //therefore reduce master/slave communication to only those IDs that 
  //require linking
  setTrackletVelocities(allDetections, allTracklets);

  std::map<ImageTime, KDTree <uint> > trackletTimeToTreeMap;
  makeTrackletTimeToTreeMap(allDetections, allTracklets, trackletTimeToTreeMap, LEAF_SIZE);
  std::cerr << timestamp() << "Worker " << rank << " made trees" << std::endl;

  //timing stuff
  int numWorkItems = 0;
  bool first = true;
  double myBFTime = 0;
  std::vector<std::vector<int> > allEndpoints;
  unsigned long long myTotalTracks = 0;

  //continually wait for tasks to be assigned
  while(true){

    MPI_Status status;

    /********************************************************
     * GET first, second and support nodes from  master node
     * as specified by filename.  Read the data from the
     * file and run brute force computation on it.
     ********************************************************/

    //receive number of TATs in this work item
    int numEndpoints;
    MPI_Recv(&numEndpoints, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    /***************************************************************************
     * numEndpoints == -1 is master node's signal to stop working
     * numEndpoints == -2 indicates a work item set assignment is complete
     * numEndpoints == 0 is an empty work item set assignment
     ***************************************************************************/
    if( numEndpoints == 0 || numEndpoints == -1 || numEndpoints == -2 ){

      double start = std::clock();

      /**************************************************
       * OPTIMAL PAGE REPLACEMENT ALGORITHM
       * 
       * Here begins the implementation of the optimal 
       * page replacement algorithm.  First we sort the
       * work items by maximum overlapping constituent
       * nodes.  Next we perform "track extraction" on 
       * those work items in the specified order, using
       * optimal page replacement to load and remove node
       * data from memory.
       **************************************************/
      
      //if no endpoints are received, skip the track extraction logic and signal
      //the master that we are finished
      if(allEndpoints.empty() || numEndpoints == 0 ){
	std::cerr << timestamp() << "There are no work item sets, continuing..." << std::endl;

	//tell the master we are finished and ready for next assignment
	int writing = 1;
	MPI_Ssend(&writing, 1, MPI_INT, 0, MPI_ANNEAL_TAG, MPI_COMM_WORLD);

	continue;
      }

      std::vector<std::vector<int> > finalEndpointOrder;

      if(DO_OPRA_1){
	finalEndpointOrder.push_back(allEndpoints.at(0));
	allEndpoints.erase(allEndpoints.begin());
      }
      else{
	for(uint i=0; i < allEndpoints.size(); ++i){
	  finalEndpointOrder.push_back(allEndpoints.at(i));
	}
	for(uint i=0; i < allEndpoints.size(); ++i){
	  allEndpoints.erase(allEndpoints.begin());
	}
      }
	
      //use allEndpoints to iterate through the work items and 
      //sort them accordingly
      if(DO_OPRA_1){
	double opra_1_start = std::clock();
	
	while(allEndpoints.size() > 0){
	  
	  //tally the number of common endpoints of each work item
	  int matches[allEndpoints.size()];

	  //initialize tally
	  for(int m=0; m < allEndpoints.size(); ++m){
	    matches[m] = 0;
	  }
	  
	  int latestWorkItemIndex = (finalEndpointOrder.size() - 1);
	  int latestWorkItemSize = finalEndpointOrder.at(latestWorkItemIndex).size();
	  
	  int bestMatch = 0;
	  int matchIndex = 0;
	  
	  for(uint i=0; i < latestWorkItemSize; i+=2){
	    
	    //compare this work item to all the remaining work items
	    for(uint k=0; k < allEndpoints.size(); ++k){
	      for(uint l=0; l < allEndpoints.at(k).size(); l+=2){
		if((finalEndpointOrder.at(latestWorkItemIndex).at(i) == allEndpoints.at(k).at(l)) &&
		   (finalEndpointOrder.at(latestWorkItemIndex).at(i+1) == allEndpoints.at(k).at(l+1))){
		  ++matches[k];
		}
	      }
	    }
	    
	    //keep track of the current best fit
	    for(uint k=0; k < allEndpoints.size(); ++k){
	      if(matches[k] > bestMatch){
		bestMatch = matches[k];
		matchIndex = k;
	      }
	    }
	  }
	
	  //update the final work item order and remove the work item from the free pool
	  finalEndpointOrder.push_back(allEndpoints.at(matchIndex));
	  allEndpoints.erase(allEndpoints.begin() + matchIndex);
	}

	double opra_1_end = getTimeElapsed(opra_1_start);
	std::cerr << timestamp() << "OPRA loop 1 took " << opra_1_end << " seconds." << std::endl;
      }

      //iterate through finalEndpointOrder vector and perform "track extraction"
      //on its work item
      //NOTE: finalEndpointOrder is a vector of work items
      double totalOpra2Time = 0;
      unsigned long long pageFaults = 0;
      for(uint i=0; i < finalEndpointOrder.size(); ++i){
	
	double opra_2_start = std::clock();

	//set up endpoint values for brute force computation
	int e1, e2;
	std::vector<treeIdNodeIdPair> supportPoints;

	for(uint j=0; j < finalEndpointOrder.at(i).size();){
      
	  //get time and ID of this endpoint
	  int treeId = finalEndpointOrder.at(i).at(j);
	  ++j;
	  int nodeId = finalEndpointOrder.at(i).at(j);
	  ++j;

	  /*************************************************************
	   * THIS IS WHERE WE WILL USE THE CACHE.  SHOULD CALL A FUNCTION
	   * THAT LOADS THIS NODE, FROM THE CACHE IF AVAILABLE, FROM
	   * DISK IF NOT.
	   *************************************************************/

	  //first iteration receives first endpoint
	  //minus 2 is because we'r incrementing the loop index inside the loop
	  if((j-2) == 0){
	    e1 = loadNodeFromCache(treeId, nodeId, cacheSize, 
				   &nodeCache[0], trackletTimeToTreeMap,
				   finalEndpointOrder, pageFaults);
	    if(e1 == -1)
	      std::cerr << timestamp() << "Firstendpoint is NULL" << std::endl;
	    }
	    //second iteration receives second endpoint
	    //minus 2 is because we'r incrementing the loop index inside the loop
	    else if((j-2) == 2){
	      e2 = loadNodeFromCache(treeId, nodeId, cacheSize, 
				     &nodeCache[0], trackletTimeToTreeMap,
				     finalEndpointOrder, pageFaults);
	      if(e2 == -1)
		std::cerr << timestamp() << "Second endpoint is NULL" << std::endl;
	    }
	    //all other iterations are support points
	    else{
	      treeIdNodeIdPair tempPair;
	      tempPair.treeId = treeId;
	      tempPair.nodeId = nodeId;
	      supportPoints.push_back(tempPair);
	    }
	}

	double opra_2_end = getTimeElapsed(opra_2_start);
	totalOpra2Time += opra_2_end;
	double workItemStart = 0, wTime = 0;

	if(e1 != -1 && e2 != -1){
	  workItemStart = std::clock();
	
	  buildTracksAddToResults(allDetections, allTracklets, searchConfig,
				  nodeCache[e1].node, nodeCache[e2].node,
				  supportPoints, 
				  myTracks, cacheSize, &nodeCache[0],
				  trackletTimeToTreeMap,
				  finalEndpointOrder, pageFaults);
	  wTime = getTimeElapsed(workItemStart);
	}
	
	//TIMING STUFF
	myBFTime += wTime;
      }
      /*
       * END OPTIMAL PAGE REPLACEMENT
       **************************************************/


      //output processing time
      double diff = getTimeElapsed(start);
      std::string thisOut("Rank ");
      std::cerr << timestamp() << thisOut << rank << " got " << numWorkItems << " work items" 
		<< ", spent " << diff << " seconds processing, " << totalOpra2Time << " seconds in OPRA 2, " 
		<< myBFTime << " seconds in brute force, and found " << myTracks.size() << " tracks with " 
		<< pageFaults << " page faults." << std::endl;

      //reset timer
      myBFTime = 0;
      
      //not returning results to master, write them to disk
      if(WRITE_RESULTS_TO_DISK){
	writeMyResults(stringify(rank), allDetections, allTracklets, myTracks, myTotalTracks);

	//tell the master we are finished and ready for next assignment
	int writing = 1;
	MPI_Ssend(&writing, 1, MPI_INT, 0, MPI_ANNEAL_TAG, MPI_COMM_WORLD);
      }
      //return results to the master
      else{
	//return all of the Tracks found, if any
	//send any track data back to the master node
	uint resultSetSize = myTracks.size();

	//inform the master how many Tracks to expect
	MPI_Ssend(&resultSetSize, 1, MPI_INT, 0, MPI_COLLECT_TAG, MPI_COMM_WORLD);
	
	//set Track data, if it exists
	if(resultSetSize > 0){
	  
	  //send data from each Track individually
	  std::set<Track>::iterator trackIter;
	  for(trackIter = myTracks.componentTracks.begin(); trackIter != myTracks.componentTracks.end(); trackIter++){
	    Track thisTrack = *trackIter;
	    
	    //send this track's componenttrackletindices
	    std::set<uint>::iterator cIt;
	    int numComponents = thisTrack.componentTrackletIndices.size();
	    int componentTrackletIndices[numComponents];
	    int count = 0;
	    for(cIt = thisTrack.componentTrackletIndices.begin();
		cIt != thisTrack.componentTrackletIndices.end();
		cIt++){
	      componentTrackletIndices[count] = *cIt; ++count;
	    }
	    //MPI_Send(&numComponents, 1, MPI_INT, 0, numComponents, MPI_COMM_WORLD);
	    //MPI_Send(&componentTrackletIndices, numComponents, MPI_INT, 0, numComponents, MPI_COMM_WORLD);
	    MPI_Send(&numComponents, 1, MPI_INT, 0, MPI_COLLECT_TAG, MPI_COMM_WORLD);
	    MPI_Send(&componentTrackletIndices, numComponents, MPI_INT, 0, MPI_COLLECT_TAG, MPI_COMM_WORLD);
	    
	    //send this track's componentdetectionindices
	    std::set<uint>::iterator dIt;
	    int numDetections = thisTrack.componentDetectionIndices.size();
	    int componentDetectionIndices[numDetections];
	    count = 0;
	    for(dIt = thisTrack.componentDetectionIndices.begin();
		dIt != thisTrack.componentDetectionIndices.end();
		dIt++){
	      componentDetectionIndices[count] = *dIt; ++count;
	    }
	    //MPI_Send(&numDetections, 1, MPI_INT, 0, numDetections, MPI_COMM_WORLD);
	    //MPI_Send(&componentDetectionIndices, numDetections, MPI_INT, 0, numDetections, MPI_COMM_WORLD);
	    MPI_Send(&numComponents, 1, MPI_INT, 0, MPI_COLLECT_TAG, MPI_COMM_WORLD);
	    MPI_Send(&componentTrackletIndices, numComponents, MPI_INT, 0, MPI_COLLECT_TAG, MPI_COMM_WORLD);
	  }
	}
	
	//tell the master that I'm quitting
	int signal = -1;
	MPI_Send(&signal, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      }

      myTotalTracks += myTracks.size();
      
      //keep going if this is just an intermediate processing block
      if( numEndpoints == -2 ){
	for(uint i=0; i < allEndpoints.size(); ++i){
	  allEndpoints.at(i).clear();
	}
	allEndpoints.clear();
	numWorkItems = 0;
	myTracks.componentTracks.clear();

	continue;
      }
      //this is the final signal, quit after this
      else if(numEndpoints == -1){
	break;
      }
      else{
	std::cerr << timestamp() << "Unexpected number of endpoints received in termination block. Aborting." << std::endl;
	break;
      }
    }
    //this is a legitimate work item -- all we do here is buffer the endpoint indices of
    //all the work items we are assigned.  Once we have been assigned all of our work items,
    //we sort them and perform "track extraction" on them in that order.
    else{
      
      //receive the node index data
      int index = 0;
      int bufferSize = numEndpoints * 2;
      int indexBuffer[bufferSize];
      MPI_Recv(&indexBuffer, bufferSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      //update the list of endpoints
      std::vector<int> temp;
      for(int i=0; i < bufferSize; ++i){
	temp.push_back(indexBuffer[i]);
      }
      
      allEndpoints.push_back(temp);
      ++numWorkItems;
    } /* end termination check */
  } /* end while(true) loop */
}






/**************************************************
 *
 **************************************************/
void distributeCurrentWorkload(std::vector<std::vector<int> > &assignment,
			       const std::vector<MopsDetection> &allDetections,
			       const std::vector<Tracklet> &allTracklets,
			       linkTrackletsConfig searchConfig,
			       std::map<ImageTime, KDTree <uint> > &trackletTimeToTreeMap,
			       /*TrackSet &results, */int currentTotalWorkItems, int currentTimeUnits,
			       std::string dir)
{
  int LEAF_SIZE = 32;
  uint numProcs = assignment.size();
  workItemFile workItemCache[WORK_ITEM_CACHE_SIZE];

  //only have one node distribution issued at a time, we do the simulated annealing
  //step if we are already waiting on threads to complete their processing. SA takes
  //a while, so may as well make that waiting useful.
  //determine which work items will be assigned to which workers
  double determineAssignmentStart = std::clock();
  determineWorkItemSets(assignment, currentTotalWorkItems, currentTimeUnits, &workItemCache[0], dir);
  double assTime = getTimeElapsed(determineAssignmentStart);
  std::cerr << timestamp() << "Master determined work item set assignment in " << assTime << " seconds." << std::endl;
 
  //find the LNS that has the minimum total time, only do this once
  /*  if(firstDistributionCycle){
      std::cerr << "Starting adaptive algorithm" << std::endl;
      double adaptiveStart = std::clock();
      int optimalLNS = adaptiveAlgorithm(allDetections, allTracklets, searchConfig, trackletTimeToTreeMap, 
      parentTable, workItemSets, LEAF_SIZE, dir);
      double adaptiveTotal = getTimeElapsed(adaptiveStart);
      std::cerr << "Adaptive algorithm took " << adaptiveTotal << " seconds." << std::endl;
      
      //send LNS to worker nodes
      for(int i=1; i <= numProcs; ++i){
      MPI_Send(&optimalLNS, 1, MPI_INT, i, i, MPI_COMM_WORLD);    
      }
      
      firstDistributionCycle = false;
      }*/
  
  //master node will distribute the work items created during the recursive tree walk
  //to the worker nodes
  distributeWorkload(assignment, &workItemCache[0], dir);
}

/***************************************************************************
 ***************************************************************************
 **   END DISTRIBUTED LINKTRACKLETS FUNCTION DEFINITIONS
 ***************************************************************************
 ***************************************************************************/




void distributedlinkTracklets(const std::vector<MopsDetection> &allDetections,
                              std::vector<Tracklet> &queryTracklets,
                              linkTrackletsConfig searchConfig, int numProcs)
{
  int LEAF_SIZE = 32;

  /*
    create a sorted list of KDtrees, each tree holding tracklets
    with unique start times (times of first detection in the tracklet).
    
    the points in the trees are in [RA, Dec, RAVelocity, DecVelocity] and the
    returned keys are indices into queryTracklets.
  */
  setTrackletVelocities(allDetections, queryTracklets);
  
  std::map<ImageTime, KDTree <uint> > trackletTimeToTreeMap;    
  makeTrackletTimeToTreeMap(allDetections, queryTracklets, trackletTimeToTreeMap, LEAF_SIZE);
  std::cerr << timestamp() << "Master made tree maps" << std::endl;

  //work item set assignment creation
  std::vector<std::vector<int> > assignment;
  std::vector<int> dummyVector;
  for(int i=0; i < (numProcs-1); ++i){
    assignment.push_back(dummyVector);
  }
  
  //create initial work item directory
  makeNextDirectory();

  //TIMING STUFF
  double treeWalkStart = std::clock();
  std::cerr << timestamp() << "Master doing linking" << std::endl;
  //std::map<uint, std::map<uint, uint> > treeIdNodeIdParentTable;
  doLinking(allDetections, queryTracklets, searchConfig, trackletTimeToTreeMap, 
	    numProcs, assignment);
  double wTime = getTimeElapsed(treeWalkStart);
  std::cerr << timestamp() << "Master did linking in " << wTime << " seconds." << std::endl;

  //catch this case here, it prevents a lot of catching throughout the code
  if(totalWorkItems > 0){
    std::cerr << timestamp() << "Distributing remaining work items." << std::endl;
    std::string dirNo = stringify(numDistributions);
    dirNo = "workItems/dir_" + dirNo;

    distributeCurrentWorkload(assignment, allDetections, queryTracklets, searchConfig, 
			      trackletTimeToTreeMap, totalWorkItems,
			      totalTimeUnits, dirNo);

    std::cerr << timestamp() << "Killing final procs" << std::endl;
    threadArgs ta;
    ta.numProcs = numProcs;
    ta.sentinel = -1;
    killProcs(ta);
  }
  //distribute any remaining work items
  else{
    std::cerr << timestamp() << "There are no more work items for this input set." << std::endl;
  }

  //find the LNS that has the minimum total time
  /*std::cerr << "Starting adaptive algorithm" << std::endl;
  double adaptiveStart = std::clock();
  int optimalLNS = adaptiveAlgorithm(allDetections, queryTracklets, searchConfig, trackletTimeToTreeMap, 
				     treeIdNodeIdParentTable, workItemSets, &workItemCache[0], LEAF_SIZE, dir);
  double adaptiveTotal = getTimeElapsed(adaptiveStart);
  std::cerr << "Adaptive algorithm took " << adaptiveTotal << " seconds." << std::endl;
  */
  //send LNS to worker nodes
  /*for(int i=1; i < numProcs; ++i){
    MPI_Send(&optimalLNS, 1, MPI_INT, i, i, MPI_COMM_WORLD);    
    }*/

}

}} // close lsst::mops

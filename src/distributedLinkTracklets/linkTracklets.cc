// -*- LSST-C++ -*-
/* jonathan myers */

#include "mpi.h"
#include <iomanip>
#include <map>
#include <gsl/gsl_multifit.h>
#include <time.h>
#include <limits.h>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <utility>
#include <string>
#include <string.h>
#include <ctime>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <signal.h>

#include "linkTracklets.h"
#include "../KDTree.h"
#include "../rmsLineFit.h"
#include "../Exceptions.h"
#include "lruCache.h"

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
  std::vector<std::pair<unsigned int, unsigned int> > endpoints;
  std::string fileName;
  int numCompatibleEndpoints;
  int numNodes;
  unsigned long long timeUnits;
}workItemFile;

typedef struct treeIdNodeId{
  unsigned int treeId;
  unsigned int nodeId;
}treeIdNodeIdPair;

unsigned int numCachedWorkItems = 0;
unsigned long int totalWorkItems = 0;
unsigned long int totalTimeUnits = 0;
unsigned int currentNumDistributions = 0;
unsigned int numDistributions = 0;
int calculatedTotal = 0;
int globalCacheSize = 0;
TrackSet toRet;
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

namespace ctExcept = collapseTracklets::exceptions;


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


class TreeNodeAndTime {
public:
  TreeNodeAndTime(){
  }
  TreeNodeAndTime(KDTree::KDTreeNode<unsigned int> * tree, ImageTime i) {
    myTree = tree;
    myTime = i;
  }
  KDTree::KDTreeNode <unsigned int> * myTree;
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
unsigned int intify(std::string s){
  
  std::istringstream i(s);
  unsigned int x;
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
			       const std::vector<Detection> &allDetections,
			       const std::vector<Tracklet> &allTracklets,
			       linkTrackletsConfig searchConfig,
			       std::map<ImageTime, KDTree::KDTree <unsigned int> > &trackletTimeToTreeMap,
			       int currentTotalWorkItems, int timeUnits, std::string dir);
  
void *killProcs(threadArgs);

int findFarthestNodeIndex(cachedNode *nodeCache, int cacheSize,
			  const std::vector<std::vector<int> > &finalEndpointOrder);

KDTree::KDTree<unsigned int> *
findTreeById(const unsigned int id, 
	     std::map<ImageTime, KDTree::KDTree<unsigned int> > &treeMap);

ImageTime findImageTimeForTree(const unsigned int id, 
			       const std::map<ImageTime, KDTree::KDTree<unsigned int> > &treeMap);

KDTree::KDTreeNode<unsigned int> *
getNodeByIDAndTime(KDTree::KDTree<unsigned int> *myTree, unsigned int nodeId);

int loadNodeFromCache(int treeId, int nodeId, int &cacheSize, cachedNode *nodeCache,
		      std::map<ImageTime, KDTree::KDTree <unsigned int> > &trackletTimeToTreeMap,
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
Detection getFirstDetectionForTracklet(const std::vector<Detection> &allDetections,
                                       const Tracklet &t) 
{
  //double start = std::clock();
    if (t.indices.size() < 1) {
        LSST_EXCEPT(collapseTracklets::exceptions::BadParameterException,
                        "linkTracklets::getFirstDetectionForTracklet called with empty tracklet.");
    }

    Detection toRet;
    bool foundOne = false;
    std::set<unsigned int>::const_iterator indexIter;
    for (indexIter = t.indices.begin(); indexIter != t.indices.end(); indexIter++) {
        Detection curDet = allDetections.at(*indexIter);
        if ((!foundOne) || 
            (toRet.getEpochMJD() > curDet.getEpochMJD())) {
            toRet = curDet;
            foundOne = true;
        }
    }

    //getFirstDetectionForTrackletTime += getTimeElapsed(start);
    return toRet;

}





void makeTrackletTimeToTreeMap(const std::vector<Detection> &allDetections,
                               const std::vector<Tracklet> &queryTracklets,
                               std::map<ImageTime, KDTree::KDTree <unsigned int> > &newMap,
			       int LEAF_SIZE)
{
    newMap.clear();
    //sort all tracklets by their first image time; make PointAndValues from
    //these so we can build a tree.
    // allTrackletPAVsMap will map from image time -> [ all tracklets starting at that image time. ]
    std::map<double, std::vector<KDTree::PointAndValue<unsigned int> > > allTrackletPAVsMap;
    allTrackletPAVsMap.clear();
    for (unsigned int i = 0; i < queryTracklets.size(); i++) {
        Detection firstDetection =  getFirstDetectionForTracklet(allDetections, queryTracklets.at(i));
        double firstDetectionTime = firstDetection.getEpochMJD();
        KDTree::PointAndValue<unsigned int> trackletPAV;
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
    unsigned int curImageID = 0;

    std::map<double, std::vector<KDTree::PointAndValue<unsigned int> > >::iterator PAVIter;
    for (PAVIter = allTrackletPAVsMap.begin(); PAVIter != allTrackletPAVsMap.end(); PAVIter++) {
        KDTree::KDTree<unsigned int> curTree(PAVIter->second, 4, LEAF_SIZE);
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
		      std::map<ImageTime, KDTree::KDTree <unsigned int> > &trackletTimeToTreeMap,
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
  KDTree::KDTree<unsigned int> *myTree = findTreeById(treeId, trackletTimeToTreeMap);
  ImageTime it = findImageTimeForTree(treeId, trackletTimeToTreeMap);
  
  //find this node in its tree
  KDTree::KDTreeNode<unsigned int> *treeNode = NULL;
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






std::set<unsigned int> allDetsInTreeNode(KDTree::KDTreeNode<unsigned int> &t,
                                            const std::vector<Detection>&allDets,
                                            const std::vector<Tracklet>&allTracklets) 
{
    std::set<unsigned int> toRet;
    std::vector<KDTree::PointAndValue<unsigned int> >::const_iterator tIter;

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
                const std::vector<Detection> &allDetections,
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




void showNumVisits(KDTree::KDTreeNode<unsigned int> *tree, unsigned int &totalNodes, unsigned int &totalVisits) 
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
void getAllDetectionsForTracklet(const std::vector<Detection> & allDetections,
                              const Tracklet &t,
                              std::vector<Detection> &detectionsForTracklet) 
{
  //double start = std::clock();
    
    detectionsForTracklet.clear();
    std::set<unsigned int>::const_iterator trackletIndexIter;

    for (trackletIndexIter = t.indices.begin();
         trackletIndexIter != t.indices.end();
         trackletIndexIter++) {
        detectionsForTracklet.push_back(allDetections.at(*trackletIndexIter));
    }

    //getAllDetectionsForTrackletTime += getTimeElapsed(start);
}



/*
 * update position and velocity given acceleration over time.
 */
void modifyWithAcceleration(double &position, double &velocity, 
                            double acceleration, double time)
{
  //double start = std::clock();
    // use good ol' displacement = vt + .5a(t^2) 
    double newPosition = position + velocity*time + .5*acceleration*(time*time);
    double newVelocity = velocity + acceleration*time;
    position = newPosition;
    velocity = newVelocity;
    //modifyWithAccelerationTime += timeSince(start);
}




inline bool positionAndVelocityRangesOverlap(double firstPositionMin, double firstPositionMax, 
                                             double firstVelocityMin, double firstVelocityMax,
                                             double secondPositionMin, double secondPositionMax,
                                             double secondVelocityMin, double secondVelocityMax)
{

  //double start = std::clock();

    bool velocityCompatible = KDTree::Common::regionsOverlap1D_unsafe(firstVelocityMin, firstVelocityMax,
                                                                      secondVelocityMin, secondVelocityMax);
    if (!velocityCompatible) {
        rejectedOnVelocity++;
        //positionAndVelocityRangesOverlapAfterAccelerationTime += timeSince(start);
        return false;
    }

    firstPositionMin = KDTree::Common::convertToStandardDegrees(firstPositionMin);
    firstPositionMax = KDTree::Common::convertToStandardDegrees(firstPositionMax);

    secondPositionMin = KDTree::Common::convertToStandardDegrees(secondPositionMin);
    secondPositionMax = KDTree::Common::convertToStandardDegrees(secondPositionMax);

    bool positionCompatible = KDTree::Common::angularRegionsOverlapSafe(firstPositionMin, firstPositionMax,
                                                                        secondPositionMin, secondPositionMax);
    if (!positionCompatible) {
        rejectedOnPosition++;
        //positionAndVelocityRangesOverlapAfterAccelerationTime += timeSince(start);
        return false;
    }
    
    wereCompatible++;
    //positionAndVelocityRangesOverlapAfterAccelerationTime += timeSince(start);
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
        throw LSST_EXCEPT(ctExcept::ProgrammerErrorException, 
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
bool areCompatible(TreeNodeAndTime  &nodeA,
                   TreeNodeAndTime  &nodeB,
                   linkTrackletsConfig searchConfig,
                   LTCache &rangeCache)
{

    double firstRAPositionMax,  firstRAPositionMin,  secondRAPositionMax,  secondRAPositionMin;
    double firstDecPositionMax, firstDecPositionMin, secondDecPositionMax, secondDecPositionMin;
    double firstRAVelocityMax,  firstRAVelocityMin,  secondRAVelocityMax,  secondRAVelocityMin;
    double firstDecVelocityMax, firstDecVelocityMin, secondDecVelocityMax, secondDecVelocityMin;

    //double start = std::clock();

    KDTree::KDTreeNode<unsigned int> * first;
    double firstTime;

    KDTree::KDTreeNode<unsigned int> * second;
    double secondTime;

    first = nodeA.myTree;
    firstTime = nodeA.myTime.getMJD();
    unsigned int firstTimeId = nodeA.myTime.getImageId();

    second = nodeB.myTree;
    secondTime = nodeB.myTime.getMJD();
    unsigned int secondTimeId = nodeB.myTime.getImageId();


    bool RACompatible = false;
    bool DecCompatible = false;

    double deltaTime = secondTime - firstTime;

    // get the bounding box for the second query region.
    secondRAPositionMax = second->getUBounds()->at(POINT_RA) + searchConfig.detectionLocationErrorThresh;
    secondRAVelocityMax = second->getUBounds()->at(POINT_RA_VELOCITY) + searchConfig.velocityErrorThresh;

    secondRAPositionMin = second->getLBounds()->at(POINT_RA) - searchConfig.detectionLocationErrorThresh;
    secondRAVelocityMin = second->getLBounds()->at(POINT_RA_VELOCITY) - searchConfig.velocityErrorThresh;
    
    secondDecPositionMax = second->getUBounds()->at(POINT_DEC) + searchConfig.detectionLocationErrorThresh;
    secondDecVelocityMax = second->getUBounds()->at(POINT_DEC_VELOCITY) + searchConfig.velocityErrorThresh;

    secondDecPositionMin = second->getLBounds()->at(POINT_DEC) - searchConfig.detectionLocationErrorThresh;
    secondDecVelocityMin = second->getLBounds()->at(POINT_DEC_VELOCITY) - searchConfig.velocityErrorThresh;

       

    // see if we have already computed the projected bounding box for the first region
    // out to this time and saved that to our cache.
    GlobalNodeIdAndProjectionId lookupKey(first->getId(), firstTimeId, secondTimeId);

    std::vector<std::vector<double> > cachedBounds;

    if (rangeCache.find(lookupKey, cachedBounds)) {
        // we don't have to recompute it -it's there already!
        cacheHits++;
        std::vector<double> uBounds = cachedBounds.at(0);
        std::vector<double> lBounds = cachedBounds.at(1);
        
        firstRAPositionMax = uBounds.at(POINT_RA);
        firstRAPositionMin = lBounds.at(POINT_RA);

        firstRAVelocityMax = uBounds.at(POINT_RA_VELOCITY);
        firstRAVelocityMin = lBounds.at(POINT_RA_VELOCITY);
        
        firstDecPositionMax = uBounds.at(POINT_DEC);
        firstDecPositionMin = lBounds.at(POINT_DEC);

        firstDecVelocityMax = uBounds.at(POINT_DEC_VELOCITY);
        firstDecVelocityMin = lBounds.at(POINT_DEC_VELOCITY);
        
        if (
            (firstRAVelocityMax < first->getUBounds()->at(POINT_RA_VELOCITY))
            ||
            (firstRAVelocityMin > first->getLBounds()->at(POINT_RA_VELOCITY))) {
                throw LSST_EXCEPT(ctExcept::ProgrammerErrorException, 
                                  "Found cached bounds at another time more restrictive than bounds at current time!");
            }
        
    }
    else {
        
        // we have not yet projected this bounding box to the given image time. 
        // do so, adding in error, and save the results.
        
        cacheMisses++;
        
        //get upper and lower bounds of node's RA position, velocity
        //modify the positional bounds using error thresholds provided by user
        /*
          TBD: we need a way to track the possible extensions to the velocity range which would
          be caused by observational error.
        */
        
        
        firstRAPositionMax = first->getUBounds()->at(POINT_RA) + searchConfig.detectionLocationErrorThresh;
        firstRAVelocityMax = first->getUBounds()->at(POINT_RA_VELOCITY) + searchConfig.velocityErrorThresh;

        firstRAPositionMin = first->getLBounds()->at(POINT_RA) - searchConfig.detectionLocationErrorThresh;
        firstRAVelocityMin = first->getLBounds()->at(POINT_RA_VELOCITY) - searchConfig.velocityErrorThresh;

        //get upper and lower bounds of node's Dec position, velocity
                
        firstDecPositionMax = first->getUBounds()->at(POINT_DEC) + searchConfig.detectionLocationErrorThresh;
        firstDecVelocityMax = first->getUBounds()->at(POINT_DEC_VELOCITY) + searchConfig.velocityErrorThresh;

        firstDecPositionMin = first->getLBounds()->at(POINT_DEC) - searchConfig.detectionLocationErrorThresh;
        firstDecVelocityMin = first->getLBounds()->at(POINT_DEC_VELOCITY) - searchConfig.velocityErrorThresh;

        double oldRAVelocityMax = firstRAVelocityMax;
        double oldRAVelocityMin = firstRAVelocityMin;
        double oldDecVelocityMax = firstDecVelocityMax;
        double oldDecVelocityMin = firstDecVelocityMin;

        if (deltaTime >= 0) {

            modifyWithAcceleration(firstRAPositionMax, firstRAVelocityMax, 
                                   searchConfig.maxRAAccel, deltaTime);
            
            modifyWithAcceleration(firstRAPositionMin, firstRAVelocityMin, 
                                   searchConfig.maxRAAccel * -1.0, deltaTime);

            modifyWithAcceleration(firstDecPositionMax, firstDecVelocityMax, 
                                   searchConfig.maxDecAccel, deltaTime);
            
            modifyWithAcceleration(firstDecPositionMin, firstDecVelocityMin, 
                                   searchConfig.maxDecAccel * -1.0, fabs(deltaTime));
            
            
        }
        else {
            // we are asking what region of RA/Dec/dRA/dDec space could
            // *reach* this point 

            extendRangeBackward(firstRAPositionMin,  firstRAPositionMax,  
                                firstRAVelocityMin,  firstRAVelocityMax,  
                                searchConfig.maxRAAccel, deltaTime);

            extendRangeBackward(firstDecPositionMin, firstDecPositionMax, 
                                firstDecVelocityMin, firstDecVelocityMax, searchConfig.maxDecAccel, deltaTime);
            
        }

        if ((oldRAVelocityMax - oldRAVelocityMin) > (firstRAVelocityMax - firstRAVelocityMin)) {
            throw LSST_EXCEPT(ctExcept::ProgrammerErrorException, "Found dec velocity range SHRUNK");
        }
        if ((oldDecVelocityMax - oldDecVelocityMin) > (firstDecVelocityMax - firstDecVelocityMin)) {
            throw LSST_EXCEPT(ctExcept::ProgrammerErrorException, "Found RA velocity range SHRUNK");
        }
        
        

        // save these newly-computed values to cache!
        std::vector<std::vector <double> > boundsForStorage(2);
        std::vector<double> uBounds(4);
        std::vector<double> lBounds(4);

        uBounds.at(POINT_RA) =  firstRAPositionMax;
        uBounds.at(POINT_DEC) = firstDecPositionMax;
        uBounds.at(POINT_RA_VELOCITY) = firstRAVelocityMax;
        uBounds.at(POINT_DEC_VELOCITY)= firstDecVelocityMax;

        lBounds.at(POINT_RA) = firstRAPositionMin;
        lBounds.at(POINT_DEC)= firstDecPositionMin;
        lBounds.at(POINT_RA_VELOCITY) = firstRAVelocityMin;
        lBounds.at(POINT_DEC_VELOCITY) = firstDecVelocityMin;

        boundsForStorage.at(0) = uBounds;
        boundsForStorage.at(1) = lBounds;

        rangeCache.insert(lookupKey, boundsForStorage);
    }
   
    
    DecCompatible = positionAndVelocityRangesOverlap(firstDecPositionMin, firstDecPositionMax,
                                                     firstDecVelocityMin, firstDecVelocityMax,
                                                     secondDecPositionMin, secondDecPositionMax,
                                                     secondDecVelocityMin, secondDecVelocityMax);
    if (!DecCompatible) {
      //areCompatibleTime += timeSince(start);
        return false;

    }

    RACompatible = positionAndVelocityRangesOverlap(firstRAPositionMin, firstRAPositionMax,
                                                    firstRAVelocityMin, firstRAVelocityMax,
                                                    secondRAPositionMin, secondRAPositionMax,
                                                    secondRAVelocityMin, secondRAVelocityMax);


    if ((RACompatible == false)) {
      //areCompatibleTime += timeSince(start);
        return false;
    }
    
    //areCompatibleTime += timeSince(start);
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
void setTrackletVelocities(const std::vector<Detection> &allDetections,
                           std::vector<Tracklet> &queryTracklets)
{
    for (unsigned int i = 0; i < queryTracklets.size(); i++) {
        Tracklet *curTracklet = &queryTracklets.at(i);
        std::vector <Detection> trackletDets;
        getAllDetectionsForTracklet(allDetections, *curTracklet, trackletDets);

        std::vector<double> RASlopeAndOffset;
        std::vector<double> DecSlopeAndOffset;
        rmsLineFit::leastSquaresSolveForRADecLinear(&trackletDets,
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
    
    //double start = std::clock();
    if (positions.size() != times.size()) {
        throw LSST_EXCEPT(ctExcept::ProgrammerErrorException,
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
    acceleration = gsl_vector_get(c,2);
    gsl_vector_free(y);
    gsl_matrix_free(X);
    gsl_vector_free(c);
    gsl_matrix_free(covariance);
    gsl_multifit_linear_free(work);
    //getBestFitVelocityAndAccelerationTime += timeSince(start);
}








void getBestFitVelocityAndAccelerationForTracklets(const std::vector<Detection> &allDetections,
                                                   const std::vector<Tracklet> &queryTracklets,
                                                   const unsigned int trackletID1,
                                                   const unsigned int trackletID2,
                                                   double & RAVelocity, double & RAAcceleration, 
                                                   double & RAPosition0,
                                                   double & DecVelocity, double & DecAcceleration, 
                                                   double & DecPosition0, 
                                                   double & time0) 
{

  //double start = std::clock();
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
        const Detection * curDetection = &(allDetections.at(*detIter));
        RAs.push_back(curDetection->getRA());
        Decs.push_back(curDetection->getDec());
        times.push_back(curDetection->getEpochMJD() - firstTime);
    }

    getBestFitVelocityAndAcceleration(RAs, times,  RAVelocity,  RAAcceleration,  RAPosition0 );
    getBestFitVelocityAndAcceleration(Decs, times, DecVelocity, DecAcceleration, DecPosition0);
    //getBestFitVelocityAndAccelerationForTrackletsTime += timeSince(start);
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
void addBestCompatibleTrackletsAndDetectionsToTrack(const std::vector<Detection> &allDetections, 
                                                    const std::vector<Tracklet> &allTracklets, 
                                                    const std::vector<unsigned int> candidateTrackletIDs, 
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
            double predRA  = RAPosition0 + RAVelocity*timeOffset + RAAcceleration*timeOffset*timeOffset;
            double predDec = DecPosition0 + DecVelocity*timeOffset + DecAcceleration*timeOffset*timeOffset;
            double distance = KDTree::Common::angularDistanceRADec_deg(detRA, detDec, predRA, predDec);
            
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
    //    addBestCompatibleTrackletsAndDetectionsToTrackTime += timeSince(start);
}







bool trackMeetsRequirements(const std::vector<Detection> & allDetections, 
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
bool endpointTrackletsAreCompatible(const std::vector<Detection> & allDetections, 
                                    const std::vector<Tracklet> &allTracklets,
                                    unsigned int trackletID1,
                                    unsigned int trackletID2,
                                    double RAVelocity, double RAAcceleration, double RAPosition0,
                                    double DecVelocity, double DecAcceleration, double DecPosition0,
                                    double time0,
                                    linkTrackletsConfig searchConfig)
{
  //double start = std::clock();
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
            double RAPred = RAPosition0 + RAVelocity*timeOffset + RAAcceleration*timeOffset*timeOffset;
            double DecPred = DecPosition0 + DecVelocity*timeOffset + DecAcceleration*timeOffset*timeOffset;
            double observedRA = allDetections.at(*detIter).getRA();
            double observedDec = allDetections.at(*detIter).getDec();
            double distanceError = 
                KDTree::Common::angularDistanceRADec_deg(RAPred, DecPred, observedRA, observedDec);
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
    //endpointTrackletsAreCompatibleTime += timeSince(start);
    return allOK;
}






/**********************************************************************
 * Matt using so each processor can write out its own results
 **********************************************************************/
void writeMyResults(std::string outFileName, 
		    const std::vector<Detection> &allDets,
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
    std::set<unsigned int>::const_iterator detIter;
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
void buildTracksAddToResults(const std::vector<Detection> &allDetections,
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
			     std::map<ImageTime, KDTree::KDTree <unsigned int> > &trackletTimeToTreeMap,
			     const std::vector<std::vector<int> > &finalEndpointOrder,
			     unsigned long long &pageFaults)
{
  //double start = std::clock();
  buildTracksAddToResultsVisits++;
  
  if ((firstEndpoint.myTree->isLeaf() == false) ||
      (secondEndpoint.myTree->isLeaf() == false)) {
    LSST_EXCEPT(ctExcept::ProgrammerErrorException, 
		"buildTracksAddToResults got non-leaf nodes, must be a bug!");
  }
  /*for (unsigned int i = 0; i < supportNodes.size(); i++) {
    if (supportNodes.at(i).myTree->isLeaf() == false) {
      LSST_EXCEPT(ctExcept::ProgrammerErrorException, 
		  "buildTracksAddToResults got non-leaf nodes, must be a bug!");            
    }
    }*/
  
  std::vector<KDTree::PointAndValue<unsigned int> >::const_iterator firstEndpointIter;
  std::vector<KDTree::PointAndValue<unsigned int> >::const_iterator secondEndpointIter;
  
  /*const std::vector<KDTree::PointAndValue<unsigned int> > *firstEndpointData = firstEndpoint.myTree->getMyData();
    const std::vector<KDTree::PointAndValue<unsigned int> > *secondEndpointData = secondEndpoint.myTree->getMyData();*/
  const std::vector<KDTree::PointAndValue<unsigned int> > *firstEndpointData = firstEndpoint.myTree->getMyData();
  const std::vector<KDTree::PointAndValue<unsigned int> > *secondEndpointData = secondEndpoint.myTree->getMyData();
  
  for (firstEndpointIter = firstEndpointData->begin();
       firstEndpointIter != firstEndpointData->end();
       firstEndpointIter++) {
    for (secondEndpointIter = secondEndpointData->begin();
	 secondEndpointIter != secondEndpointData->end();
	 secondEndpointIter++) {

      //std::vector<TreeNodeAndTime>::const_iterator supportNodeIter;
      //std::vector<KDTree::PointAndValue <unsigned int> >::const_iterator supportPointIter;
      std::vector<treeIdNodeIdPair>::const_iterator supportNodeIter;
      std::vector<KDTree::PointAndValue <unsigned int> >::const_iterator supportPointIter;
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
	//for (supportNodeIter = supportNodes.begin(); supportNodeIter != supportNodes.end();
	//supportNodeIter++) {
	for (supportNodeIter = supportNodes.begin(); supportNodeIter != supportNodes.end();
	     supportNodeIter++) {
	  treeIdNodeIdPair tini = *supportNodeIter;
	  int i = loadNodeFromCache(tini.treeId, tini.nodeId, cacheSize, 
				    nodeCache, trackletTimeToTreeMap,
				    finalEndpointOrder, pageFaults);
	  
	  const std::vector<KDTree::PointAndValue <unsigned int> > * curSupportNodeData;
	  //if(!nodeCache[i].node->myTree->isLeaf()){
	  if(!nodeCache[i].node.myTree->isLeaf()){
	    throw LSST_EXCEPT(ctExcept::BadParameterException,
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
  unsigned int count = 0;
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







double nodeWidth(KDTree::KDTreeNode<unsigned int> *node)
{
  double width = 1;
  for (unsigned int i = 0; i < 4; i++) {
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
KDTree::KDTree<unsigned int> *
findTreeById(const unsigned int id, 
	     std::map<ImageTime, KDTree::KDTree<unsigned int> > &treeMap)
{
  KDTree::KDTree<unsigned int> *retVal = NULL;

  std::map<ImageTime, KDTree::KDTree<unsigned int> >::iterator mapIter;
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
ImageTime findImageTimeForTree(const unsigned int id, 
			       const std::map<ImageTime, KDTree::KDTree<unsigned int> > &treeMap)
{

  std::map<ImageTime, KDTree::KDTree<unsigned int> >::const_iterator mapIter;
  for(mapIter = treeMap.begin(); mapIter != treeMap.end(); mapIter++){
    ImageTime it = (*mapIter).first;
    if(it.getImageId() == id){
      return it;
    }
  }

  ImageTime retVal;
  return retVal;
}



KDTree::KDTreeNode<unsigned int> *
getNodeByIDAndTime(KDTree::KDTree<unsigned int> *myTree, unsigned int nodeId)
{
  KDTree::KDTreeNode<unsigned int> *retVal;

  //search tree for nodeId
  std::vector<KDTree::KDTreeNode<unsigned int> *> nodeList;
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
void doLinkingRecurse2(const std::vector<Detection> &allDetections,
                       const std::vector<Tracklet> &allTracklets,
                       linkTrackletsConfig searchConfig,
                       TreeNodeAndTime &firstEndpoint,
                       TreeNodeAndTime &secondEndpoint,
                       std::vector<TreeNodeAndTime> &supportNodes,
                       //TrackSet & results,
                       int iterationsTillSplit,
                       LTCache &rangeCache, 
		       int numProcs, std::map<ImageTime, KDTree::KDTree <unsigned int> > &trackletTimeToTreeMap,
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
		    unsigned int treeId = supportNodeIter->myTime.getImageId();
		    unsigned int parentId = supportNodeIter->myTree->getId();
		    unsigned int childId = newTAT.myTree->getId();
		    parentTable[treeId][childId] = parentId;
		  */
		  newSupportNodes.push_back(newTAT);
		  uniqueSupportMJDs.insert(supportNodeIter->myTime.getMJD());
		}
		if (supportNodeIter->myTree->hasRightChild()) {
		  
		  TreeNodeAndTime newTAT(supportNodeIter->myTree->getRightChild(), supportNodeIter->myTime);
		  /*
		    unsigned int treeId = firstEndpoint.myTime.getImageId();
		    unsigned int parentId = firstEndpoint.myTree->getId();
		    unsigned int childId = newTAT.myTree->getId();
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
	    
	    std::vector<KDTree::PointAndValue<unsigned int> >::const_iterator firstEndpointIter;
	    std::vector<KDTree::PointAndValue<unsigned int> >::const_iterator secondEndpointIter;
	    
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
	    std::vector<std::pair<unsigned int, unsigned int> > _endpoints;
	    _endpoints.clear();
	    
	    //first item is first endpoint
	    unsigned int treeId = firstEndpoint.myTime.getImageId();
	    unsigned int nodeId = firstEndpoint.myTree->getId();
	    _endpoints.push_back(std::make_pair(treeId, nodeId));
	    
	    //second item is second endpoint
	    treeId = secondEndpoint.myTime.getImageId();
	    nodeId = secondEndpoint.myTree->getId();
	    _endpoints.push_back(std::make_pair(treeId, nodeId));
	    
	    //add all endpoints
	    for(unsigned int m=0; m < newSupportNodes.size(); ++m){
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
	    
	    for (unsigned int i = 0; i < myItem.endpoints.size(); i++) {
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
	      for(unsigned int i=0; i < assignment.size(); ++i){
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
                if ( (KDTree::Common::areEqual(firstEndpointWidth, -1)) &&
                     (KDTree::Common::areEqual(secondEndpointWidth, -1)) ) {
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
		    throw LSST_EXCEPT(ctExcept::ProgrammerErrorException, "Recursing in a leaf node (first endpoint), must be a bug!");
		  }
		  
		  if (firstEndpoint.myTree->hasLeftChild())
                    {
		      TreeNodeAndTime newTAT(firstEndpoint.myTree->getLeftChild(), firstEndpoint.myTime);
			      
		      /**************************************************
		       *MATT UPDATING FOR ADAPTIVE ALGORITHM
		       **************************************************/
		      /*unsigned int treeId = firstEndpoint.myTime.getImageId();
			unsigned int parentId = firstEndpoint.myTree->getId();
			unsigned int childId = newTAT.myTree->getId();
			parentTable[treeId][childId] = parentId;
		      */
		      /* 
		       * END
		       **************************************************/
		      doLinkingRecurse2(allDetections, allTracklets, searchConfig,
					newTAT, secondEndpoint,
					newSupportNodes,
					iterationsTillSplit, rangeCache, 
					numProcs, trackletTimeToTreeMap, assignment);
                    }
		  
		  if (firstEndpoint.myTree->hasRightChild())
                    {
		      TreeNodeAndTime newTAT(firstEndpoint.myTree->getRightChild(), firstEndpoint.myTime);
			      
		      /**************************************************
		       *MATT UPDATING FOR ADAPTIVE ALGORITHM
		       **************************************************/
		      /*unsigned int treeId = firstEndpoint.myTime.getImageId();
		      unsigned int parentId = firstEndpoint.myTree->getId();
		      unsigned int childId = newTAT.myTree->getId();
		      parentTable[treeId][childId] = parentId;
		      */
		      /* 
		       * END
		       **************************************************/
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
		    throw LSST_EXCEPT(ctExcept::ProgrammerErrorException, "Recursing in a leaf node (second endpoint), must be a bug!");
		  }
		  
		  if (secondEndpoint.myTree->hasLeftChild())
                    {
		      TreeNodeAndTime newTAT(secondEndpoint.myTree->getLeftChild(), secondEndpoint.myTime);
			      
		      /**************************************************
		       *MATT UPDATING FOR ADAPTIVE ALGORITHM
		       **************************************************/
		      /*unsigned int treeId = secondEndpoint.myTime.getImageId();
			unsigned int parentId = secondEndpoint.myTree->getId();
			unsigned int childId = newTAT.myTree->getId();
			parentTable[treeId][childId] = parentId;
		      */
		      /* 
		       * END
		       **************************************************/
		      doLinkingRecurse2(allDetections, allTracklets, searchConfig,
					firstEndpoint,newTAT,
					newSupportNodes,
					iterationsTillSplit, rangeCache, 
					numProcs, trackletTimeToTreeMap, assignment);
                    }
		  
		  if (secondEndpoint.myTree->hasRightChild())
                    {
		      TreeNodeAndTime newTAT(secondEndpoint.myTree->getRightChild(), secondEndpoint.myTime);
		      
		      /**************************************************
		       *MATT UPDATING FOR ADAPTIVE ALGORITHM
		       **************************************************/
		      /*unsigned int treeId = secondEndpoint.myTime.getImageId();
			unsigned int parentId = secondEndpoint.myTree->getId();
			unsigned int childId = newTAT.myTree->getId();
			parentTable[treeId][childId] = parentId;
		      */
		      /* 
		       * END
		       **************************************************/
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



void doLinking(const std::vector<Detection> &allDetections,
               std::vector<Tracklet> &allTracklets,
               linkTrackletsConfig searchConfig,
               std::map<ImageTime, KDTree::KDTree <unsigned int> > &trackletTimeToTreeMap,
	       int numProcs, std::vector<std::vector<int> > &assignment)

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
  double limitedRunFirstEndpoint = 53738.190801;//53738.217607;///53738.266849;/////49616.273436000003;//   49616.249649999998 ;
  double limitedRunSecondEndpoint = 53747.220188;//53747.208854;//53747.241076;//49623.023787999999 ;  //49623.023787999999 ;
  
  double iterationTime = std::clock();
  int threshold = 1000;
  
  if (DEBUG) {
    std:: cout << "all MJDs: ";
    std::map<ImageTime, KDTree::KDTree<unsigned int> >::const_iterator mapIter;
    for (mapIter = trackletTimeToTreeMap.begin();
	 mapIter != trackletTimeToTreeMap.end();
	 mapIter++) {
      std::cout << std::setprecision (10) << mapIter->first.getMJD() << " ";
    }
    std::cout << std::endl;
  }
  
  std::map<ImageTime, KDTree::KDTree<unsigned int> >::const_iterator firstEndpointIter;
  for (firstEndpointIter = trackletTimeToTreeMap.begin(); 
       firstEndpointIter != trackletTimeToTreeMap.end(); 
       firstEndpointIter++)
    {
      std::map<ImageTime, KDTree::KDTree<unsigned int> >::const_iterator secondEndpointIter;
      std::map<ImageTime, KDTree::KDTree<unsigned int> >::const_iterator afterFirstIter = firstEndpointIter;
      afterFirstIter++;
      
      if ((!limitedRun) || (KDTree::Common::areEqual(firstEndpointIter->first.getMJD(), limitedRunFirstEndpoint))) {
	
	for (secondEndpointIter = afterFirstIter; 
	     secondEndpointIter != trackletTimeToTreeMap.end(); 
	     secondEndpointIter++)
	  {
	    /* if there is sufficient time between the first and second nodes, then try
	       doing the recursive linking using intermediate times as support nodes.
	    */
	    int numAdded = 0;
	    if ((!limitedRun) || 
		(KDTree::Common::areEqual(secondEndpointIter->first.getMJD(), limitedRunSecondEndpoint))) {
                    
	      
	      if (secondEndpointIter->first.getMJD() - firstEndpointIter->first.getMJD() 
		  >= searchConfig.minEndpointTimeSeparation) {                
		
		if (DEBUG) {
		  time_t rawtime;
                  
		  if (timeSince(iterationTime) > 5) {
		    //debugPrintTimingInfo(results);
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
		std::map<ImageTime, KDTree::KDTree<unsigned int> >::const_iterator supportPointIter;
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
		       searchConfig.minSupportToEndpointTimeSeparation)
		      && numAdded < threshold) {
		    
		    TreeNodeAndTime tmpTAT(supportPointIter->second.getRootNode(),
					   supportPointIter->first);
		    supportPoints.push_back(tmpTAT);
		    ++numAdded;
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
		//LTCache rangeCache(100000000);
		LTCache rangeCache(0);
		doLinkingRecurse2(allDetections, allTracklets, searchConfig,
				  firstEndpoint, secondEndpoint,
				  supportPoints,
				  2, rangeCache, 
				  numProcs, trackletTimeToTreeMap, assignment);
	      }
	    }
	  }
      }
    }
}







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
      
      std::set<unsigned int> tSet;
      for(int k=0; k < numTracklets; ++k){
	tSet.insert(trackletIndices[k]);
      }
      
      //receive detection indices
      int numDetections;
      MPI_Recv(&numDetections, 1, MPI_INT, _from, MPI_COLLECT_TAG, MPI_COMM_WORLD, &status);    
      
      int detectionIndices[numDetections];
      MPI_Recv(&detectionIndices, numDetections, MPI_INT, _from, MPI_COLLECT_TAG, MPI_COMM_WORLD, &status);
      
      std::set<unsigned int> dSet;
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
  unsigned int ti;
  unsigned int ni;
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
    for(unsigned int i=0; i < numCachedWorkItems; ++i){
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
  unsigned int numWorkers = assignment.size();
  for(unsigned int i=0; i < numWorkers; ++i){
    mixedSet.push_back(assignment.at(i));
  }

  //shuffle all the work items assigned to each worker
  std::random_shuffle(mixedSet.begin(), mixedSet.end());
  for(unsigned int i=0; i < mixedSet.size(); ++i){
    std::random_shuffle(mixedSet.at(i).begin(), mixedSet.at(i).end());
  }

  unsigned long long timePerWorker = (currentTimeUnits / numWorkers);
  int nextWorkItem = 0;
  unsigned int workItemNodeSize = mixedSet.at(nextWorkItem).size();
  unsigned int numAssignedWorkItems = 0;
  unsigned int workItemNodeIndex = 0;
  unsigned long long assignedTime = 0;

  //each worker gets their share of work items
  for(unsigned int i=0; i < numWorkers; ++i){
    
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
  unsigned int numRemaining = totalWorkItems - numAssignedWorkItems;

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

    for(unsigned int i = n; i < mixedSet.size(); ++i){
      for(unsigned int j = remainder; j < mixedSet.at(i).size(); ++j){
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
  unsigned int L = 1000;

  //this is the max cost from this work item set assignment
  unsigned long long maxCost = 0;

  for(unsigned int i=0; i < assignment.size(); ++i){
    
    unsigned long long cost = 0;

    //look through all the nodes of the adjacent work item and find number of shared nodes
    for(unsigned int j=0; j < assignment.at(i).size(); ++j){
      cost += executionModelCost(assignment.at(i).at(j), workItemCache, dir);
      
      workItemFile oneTemp = readWorkItemFileFromDisk(assignment.at(i).at(j), workItemCache, dir);
      unsigned int matches = 0;
      unsigned int numNodes = oneTemp.numNodes;;
      
      // memory (IO) model cost (this only goes here because the IO
      // cost varies at the WORK ITEM SET ASSIGNMENT level, not the work
      // item level)
      if(j == 0){
	cost+=(numNodes * L);
      }
      else{

	//think: should this be (j+1)?
	workItemFile twoTemp = readWorkItemFileFromDisk(assignment.at(i).at((j-1)), workItemCache, dir);

	for(unsigned int k=0; k <  oneTemp.endpoints.size(); ++k){
	  for(unsigned int l=0; l < twoTemp.endpoints.size(); ++l){
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
bool updateCurrent(unsigned int cost, unsigned int bestCost, unsigned int T)
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
  for(unsigned int i=0; i < assignment.size(); ++i){
    currentAss.push_back(assignment.at(i));
  }

  //number of iterations to consider different work item set 
  //assignment permutations
  int maxIterations = 10, throttle = 500;
  int iteration = 0;
  unsigned int T = 10;

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
      for(unsigned int i=0; i < assignment.size(); ++i){
	assignment.at(i).clear();
      }
      assignment.clear();
      
      //replace current best with the new best
      for(unsigned int i=0; i < trialAss.size(); ++i){
	assignment.push_back(trialAss.at(i));
      }
    }

    //check if we should update the current set assignment
    if(updateCurrent(trialCost, bestCost, T)){

      //clear the current set assignment
      for(unsigned int i=0; i < currentAss.size(); ++i){
	currentAss.at(i).clear();
      }
      currentAss.clear();
      
      //replace current set assignment with trial set assignment
      for(unsigned int i=0; i < trialAss.size(); ++i){
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
  unsigned int numProcs = workItemSets.size();

  //std::cerr << timestamp() << "Master distributing " << totalWorkItems << " work items to " << numProcs << " workers."  <<  std::endl;

  long int workUnitsList[numProcs];
  for(int i=0; i<numProcs; ++i){
    workUnitsList[i] = 0;
  }

  for(unsigned int i=0; i < workItemSets.size(); ++i){

    //corner case where a worker wasn't assigned any work items, inform him that he 
    //should expect zero items
    if(workItemSets.at(i).size() == 0){
      int numEndpoints = 0;
      MPI_Send(&numEndpoints, 1, MPI_INT, (i+1), i, MPI_COMM_WORLD);
    }
    else{
      for(unsigned int j=0; j < workItemSets.at(i).size(); ++j){
	
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
  
  for(unsigned int i=1; i < ta.numProcs; ++i){
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
    for(unsigned int j=0; j< finalEndpointOrder.size(); ++j){
      
      //go through this work item's nodes
      for(unsigned int k=0; k < finalEndpointOrder.at(j).size() && !found; ++k){

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
		 const std::vector<Detection> &allDetections,
		 std::vector<Tracklet> &allTracklets, 
		 linkTrackletsConfig searchConfig)
{
  //this is the set of tracks this function returns
  TrackSet myTracks;

  //this is THE node cache, declare and initialize it
  cachedNode nodeCache[CACHE_SIZE];
  for(unsigned int i=0; i < CACHE_SIZE; ++i){
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

  std::map<ImageTime, KDTree::KDTree <unsigned int> > trackletTimeToTreeMap;
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
	for(unsigned int i=0; i < allEndpoints.size(); ++i){
	  finalEndpointOrder.push_back(allEndpoints.at(i));
	}
	for(unsigned int i=0; i < allEndpoints.size(); ++i){
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
	  
	  for(unsigned int i=0; i < latestWorkItemSize; i+=2){
	    
	    //compare this work item to all the remaining work items
	    for(unsigned int k=0; k < allEndpoints.size(); ++k){
	      for(unsigned int l=0; l < allEndpoints.at(k).size(); l+=2){
		if((finalEndpointOrder.at(latestWorkItemIndex).at(i) == allEndpoints.at(k).at(l)) &&
		   (finalEndpointOrder.at(latestWorkItemIndex).at(i+1) == allEndpoints.at(k).at(l+1))){
		  ++matches[k];
		}
	      }
	    }
	    
	    //keep track of the current best fit
	    for(unsigned int k=0; k < allEndpoints.size(); ++k){
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
      for(unsigned int i=0; i < finalEndpointOrder.size(); ++i){
	
	double opra_2_start = std::clock();

	//set up endpoint values for brute force computation
	int e1, e2;
	std::vector<treeIdNodeIdPair> supportPoints;

	for(unsigned int j=0; j < finalEndpointOrder.at(i).size();){
      
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
	unsigned int resultSetSize = myTracks.size();

	//inform the master how many Tracks to expect
	MPI_Ssend(&resultSetSize, 1, MPI_INT, 0, MPI_COLLECT_TAG, MPI_COMM_WORLD);
	
	//set Track data, if it exists
	if(resultSetSize > 0){
	  
	  //send data from each Track individually
	  std::set<Track>::iterator trackIter;
	  for(trackIter = myTracks.componentTracks.begin(); trackIter != myTracks.componentTracks.end(); trackIter++){
	    Track thisTrack = *trackIter;
	    
	    //send this track's componenttrackletindices
	    std::set<unsigned int>::iterator cIt;
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
	    std::set<unsigned int>::iterator dIt;
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
	for(unsigned int i=0; i < allEndpoints.size(); ++i){
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



/*********************************************************************************
 * 
 * ADAPTIVE LNS ALGORITHM 
 * 
 *********************************************************************************/
int adaptiveAlgorithm(const std::vector<Detection> &allDetections,
		      std::vector<Tracklet> &allTracklets,
		      linkTrackletsConfig searchConfig,
		      std::map<ImageTime, KDTree::KDTree <unsigned int> > &trackletTimeToTreeMap,
		      std::map<unsigned int, std::map<unsigned int, unsigned int> > &parentTable,
		      std::vector<std::vector<int> > *assignment, workItemFile *workItemCache,
		      int LEAF_SIZE, std::string dir)
{
  std::cerr << "Entering adaptive algorithm" << std::endl;

  //this was a #define but i'm trying to minimize globals that aren't needed by both
  //master and workers
  int MAX_ADAPTIVE_LNS = 64;

  //create a map of LNS to work item vectors
  std::map<unsigned int, std::vector<int>*> workItemMap;
  int currentLNS = LEAF_SIZE;
  workItemMap[currentLNS] = &assignment->at(0);

  //keep track of timeUnits at each LNS
  std::map<int, long unsigned int> LNSTimeMap;

  //this multimap stores endpoint pairs that have already been visited
  std::vector<std::pair<treeIdNodeIdPair, treeIdNodeIdPair> > endpointsVisited;
  endpointsVisited.clear();

  //endpoint nodes
  TreeNodeAndTime *firstEndpoint;
  TreeNodeAndTime *secondEndpoint;
  std::vector<TreeNodeAndTime*> newSupportNodes;

  //variables for work item creation
  double RAVelocity, DecVelocity, RAAcceleration, DecAcceleration;
  double RAPosition0, DecPosition0;
  double time0;
  int numCompatible = 0;

  //perform time aggregation on all LNSes
  while(currentLNS <= MAX_ADAPTIVE_LNS){

    std::cerr << "Adaptive algorithm evaluating LNS of " << currentLNS << std::endl;
    if(workItemMap.find(currentLNS) != workItemMap.end()){
      
      unsigned int count = 0;
      //gather work item time for each work item of this LNS
      for(unsigned int i=0; i < workItemMap[currentLNS]->size(); ++i){
	
	//load this work item
	workItemFile adaptiveItem = readWorkItemFileFromDisk(workItemMap[currentLNS]->at(i), workItemCache, dir);

	//nodeId and treeId of the first endpoint
	//unsigned int treeId = workItemMap[currentLNS]->at(i).endpoints.at(0).first;
	//int nodeId = workItemMap[currentLNS]->at(i).endpoints.at(0).second;
	unsigned int treeId = adaptiveItem.endpoints.at(0).first;
	int nodeId = adaptiveItem.endpoints.at(0).second;
	
	//Get the tree and ImageTime for this image (tree) ID
	KDTree::KDTree<unsigned int> *myTree = findTreeById(treeId, trackletTimeToTreeMap);
	ImageTime it = findImageTimeForTree(treeId, trackletTimeToTreeMap);
	
	//find this node in its tree
	KDTree::KDTreeNode<unsigned int> *treeNode = NULL;
	if( myTree != NULL ){
	  treeNode = getNodeByIDAndTime(myTree, nodeId);
	}
	
	//ensure proper node was found
	if( treeNode != NULL ){
	  firstEndpoint = new TreeNodeAndTime(treeNode, it);
	}	
	
	//nodeId and treeId of second endpoint
	//treeId = workItemMap[currentLNS]->at(i).endpoints.at(1).first;
	//nodeId = workItemMap[currentLNS]->at(i).endpoints.at(1).second;
	treeId = adaptiveItem.endpoints.at(1).first;
	nodeId = adaptiveItem.endpoints.at(1).second;
	
	//Get the tree and ImageTime for this image (tree) ID
	myTree = findTreeById(treeId, trackletTimeToTreeMap);
	it = findImageTimeForTree(treeId, trackletTimeToTreeMap);
	
	//find this node in its tree
	treeNode = NULL;
	if( myTree != NULL ){
	  treeNode = getNodeByIDAndTime(myTree, nodeId);
	}
	
	//ensure proper node was found
	if( treeNode != NULL ){
	  secondEndpoint = new TreeNodeAndTime(treeNode, it);
	}	
	
	//gather up support points
	//for(unsigned int j=2; j < workItemMap[currentLNS]->at(i).endpoints.size(); ++j){
	for(unsigned int j=2; j < adaptiveItem.endpoints.size(); ++j){
	  
	  //treeId = workItemMap[currentLNS]->at(i).endpoints.at(j).first;
	  //nodeId = workItemMap[currentLNS]->at(i).endpoints.at(j).second;
	  treeId = adaptiveItem.endpoints.at(j).first;
	  nodeId = adaptiveItem.endpoints.at(j).second;
	  
	  //Get the tree and ImageTime for this image (tree) ID
	  myTree = findTreeById(treeId, trackletTimeToTreeMap);
	  it = findImageTimeForTree(treeId, trackletTimeToTreeMap);
	  
	  //find this node in its tree
	  treeNode = NULL;
	  if( myTree != NULL ){
	    treeNode = getNodeByIDAndTime(myTree, nodeId);
	  }
	  
	  //ensure proper node was found
	  if( treeNode != NULL ){
	    TreeNodeAndTime *TAT = new TreeNodeAndTime(treeNode, it);
	    newSupportNodes.push_back(TAT);
	  }	
	}
	
	std::vector<KDTree::PointAndValue<unsigned int> >::const_iterator firstEndpointIter;
	std::vector<KDTree::PointAndValue<unsigned int> >::const_iterator secondEndpointIter;
	const std::vector<KDTree::PointAndValue<unsigned int> > *firstEndpointData = firstEndpoint->myTree->getMyData();
	const std::vector<KDTree::PointAndValue<unsigned int> > *secondEndpointData = secondEndpoint->myTree->getMyData();
	
	for (firstEndpointIter = firstEndpointData->begin();
	     firstEndpointIter != firstEndpointData->end();
	     firstEndpointIter++) {
	  
	  for (secondEndpointIter = secondEndpointData->begin();
	       secondEndpointIter != secondEndpointData->end();
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
	
	//calculate the amount of time this work item requires and update LNSTimeMap
	//workItemMap[currentLNS]->at(i).timeUnits = ((2 + newSupportNodes.size()) * numCompatible);
	adaptiveItem.timeUnits = ((2 + newSupportNodes.size()) * numCompatible);
	
	/*if(LNSTimeMap.find(currentLNS) == LNSTimeMap.end()){
	  LNSTimeMap[currentLNS] = workItemMap[currentLNS]->at(i).timeUnits;
	}
	else{
	  LNSTimeMap[currentLNS] += workItemMap[currentLNS]->at(i).timeUnits;
	  }*/

	if(LNSTimeMap.find(currentLNS) == LNSTimeMap.end()){
	  LNSTimeMap[currentLNS] = adaptiveItem.timeUnits;
	}
	else{
	  LNSTimeMap[currentLNS] += adaptiveItem.timeUnits;
	}
	
	/******************************************************
	 * create work item for this work item's parent points
	 ******************************************************/
	//these endpoints are those for the work item of this work item's parent
	std::vector<std::pair<unsigned int, unsigned int> > endpoints;
	endpoints.clear();
	
	//first item is first endpoint
	treeId = firstEndpoint->myTime.getImageId();
	nodeId = firstEndpoint->myTree->getId();
	unsigned int parentNodeId = parentTable[treeId][nodeId];
	endpoints.push_back(std::make_pair(treeId, parentNodeId));
	
	//second item is second endpoint
	treeId = secondEndpoint->myTime.getImageId();
	nodeId = secondEndpoint->myTree->getId();
	parentNodeId = parentTable[treeId][nodeId];
	endpoints.push_back(std::make_pair(treeId, parentNodeId));
	
	//add all endpoints
	for(unsigned int m=0; m < newSupportNodes.size(); ++m){
	  treeId = newSupportNodes.at(m)->myTime.getImageId();
	  nodeId = newSupportNodes.at(m)->myTree->getId();
	  parentNodeId = parentTable[treeId][nodeId];
	  endpoints.push_back(std::make_pair(treeId, parentNodeId));
	}
	
	//create workItemFile struct to go on list
	workItemFile myItem;
	myItem.endpoints = endpoints;
	myItem.numNodes = 2 + newSupportNodes.size();
	myItem.numCompatibleEndpoints = 0;                 //numCompatible  -- we don't have this value yet
	myItem.timeUnits = myItem.numNodes * myItem.numCompatibleEndpoints;  //we don't have this value yet
	
	//clear this vector -- it gets big and is workItem-specific
	newSupportNodes.clear();
	
	//add work item to list at corresponding LNS slot
	int parentLNS = currentLNS * 2;
	if(workItemMap.find(parentLNS) == workItemMap.end()){
	  //workItemMap[parentLNS] = new std::vector<workItemFile>();
	  workItemMap[parentLNS] = new std::vector<int>();
	  workItemMap[parentLNS]->clear();
	}
	
	//std::cerr << "checking to see if this endpoint pair exists" << std::endl;
	
	//check to see if this work item already exists, don't add duplicate
	bool exists = false;
	for(unsigned int k=0; k < endpointsVisited.size(); ++k){
	  if((endpointsVisited.at(k).first.treeId == 
	      myItem.endpoints.at(0).first) &&
	     (endpointsVisited.at(k).first.nodeId ==
	      myItem.endpoints.at(0).second)){
	    
	    //only check second endpoint if first was found
	    if((endpointsVisited.at(k).second.treeId ==
		myItem.endpoints.at(1).first) &&
	       (endpointsVisited.at(k).second.nodeId ==
		myItem.endpoints.at(1).second)){
	      exists = true;
	      break;
	    }
	  }
	}
	
	//add this work item to the list if it hasn't already been added
	if(!exists){	  
	  //std::cerr << "doesn't exist, adding to endpointsVisited" << std::endl;
	  
	  //update the list of workitem endpoints visited
	  treeIdNodeIdPair p1;
	  p1.treeId = myItem.endpoints.at(0).first;
	  p1.nodeId = myItem.endpoints.at(0).second;

	  treeIdNodeIdPair p2;
	  p2.treeId = myItem.endpoints.at(1).first;
	  p2.nodeId = myItem.endpoints.at(1).second;

	  endpointsVisited.push_back(std::make_pair(p1, p2));
	  
	  //add this new work item  to thelist associated with the parents
	  //of the  current  work item
	  //workItemMap[parentLNS]->push_back(myItem);
	  workItemMap[parentLNS]->push_back(parentTable[treeId][nodeId]);
	}
	//std::cerr << "work item loop finished iteration " << i << std::endl;
	count = i;
      } /* end work item loop */ 
      std::cerr << "LNS " << currentLNS << " had " << count << " work items!?!" << std::endl;
    } /* LNS exists check */
    
    //this list is unique to each LNS
    endpointsVisited.clear();

    //update the LNS to the next higher power of 2
    currentLNS *= 2;
  } /* end adapative algorithm loop */

  std::cerr << "Finished adaptive loop" << std::endl;

  //find the LNS with lowest total work units and return it
  long unsigned int minWorkUnits = ULONG_MAX;
  int optimalLNS = LEAF_SIZE;
  bool first = true;

  std::map<int, long unsigned int>::const_iterator LNSiter;
  
  for(LNSiter = LNSTimeMap.begin(); LNSiter != LNSTimeMap.begin(); LNSiter++){
    if(((*LNSiter).second < minWorkUnits) || first){
      minWorkUnits = (*LNSiter).second;
      optimalLNS = (*LNSiter).first;
      std::cerr << "updating minWorkUnits to " << minWorkUnits << std::endl;
      std::cerr << "optimalLNS is now " << optimalLNS << std::endl;

      if(first){
	first = false;
      }
    }
  }

  //print adaptive algorithm results
  std::cerr << "LNS Total Work times" << std::endl;
  std::cerr << "----------------------------------------" << std::endl;
  std::map<int, unsigned long int>::iterator it;
  for (it = LNSTimeMap.begin() ; it != LNSTimeMap.end(); it++ )
    std::cerr << (*it).first << " => " << (*it).second << std::endl;

  return optimalLNS;
} 



/**************************************************
 *
 **************************************************/
void distributeCurrentWorkload(std::vector<std::vector<int> > &assignment,
			       const std::vector<Detection> &allDetections,
			       const std::vector<Tracklet> &allTracklets,
			       linkTrackletsConfig searchConfig,
			       std::map<ImageTime, KDTree::KDTree <unsigned int> > &trackletTimeToTreeMap,
			       /*TrackSet &results, */int currentTotalWorkItems, int currentTimeUnits,
			       std::string dir)
{
  int LEAF_SIZE = 32;
  unsigned int numProcs = assignment.size();
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




TrackSet linkTracklets(const std::vector<Detection> &allDetections,
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
  
  std::map<ImageTime, KDTree::KDTree <unsigned int> > trackletTimeToTreeMap;    
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
  //std::map<unsigned int, std::map<unsigned int, unsigned int> > treeIdNodeIdParentTable;
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

  return (toRet);
}


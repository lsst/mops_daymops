// -*- LSST-C++ -*-
/* jonathan myers and Matthew Cleveland */

#include "mpi.h"
// time headers needed for benchmarking performance
#include <sys/stat.h>
#include <sys/types.h>
#include <ctime>
#include <iomanip>
#include <map>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <limits>
#include <limits.h>

#include "lsst/mops/rmsLineFit.h"
#include "lsst/mops/daymops/linkTracklets/linkTracklets.h"
#include "lsst/mops/Exceptions.h"
#include "lsst/mops/KDTree.h"
#include "lsst/mops/daymops/linkTracklets/TrackletTree.h"



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
#define LEAF_SIZE 32


#define uint unsigned int 

/* open namespaces */
namespace lsst {
    namespace mops {



/**************************************************
 **************************************************
 **  DISTRIBUTED VARIABLE DEFINITIONS
 **************************************************
 **************************************************/
unsigned int numGenerations  = 0;
unsigned long int totalWorkItems  = 0;
unsigned long int totalTimeUnits  = 0;
unsigned long int globalCacheSize = 0;
TrackSet        * toRet;
pthread_t         thread;

//used by workers to track their total cache statistics
unsigned long int globalCacheHits = 0, globalCacheMisses = 0;

#define MAX_WORK_ITEMS        10000     /* number of work items per distribution batch */

#define MPI_NODE_COMM_TAG     9999   /* MPI tage for node communication */
#define MPI_COLLECT_TAG       6969   /* MPI tag for track collection */
#define MPI_ANNEAL_TAG        4242   /* MPI tag for annealing trigger */

#define CACHE_SIZE            64     /* size of node cache on worker processors */
#define MAX_GENERATIONS       10000  /* number of work item generations to complete */
#define WRITE_RESULTS_TO_DISK 1      /* workers write tracks to disk or return them */

#define SIMULATE_DISK_ACCESS  0      /* read random data from disk to simulate IO */
#define DO_LRU                0      /* perform least recently used cache replacment algorithm */
#define DO_ANNEALING          1      /* perform simulated annealing or not */
#define DO_OPRA               1      /* perform optimal page replacement algorithm for node cache loading */
#define DO_OPRA_PREP          0      /* perform the optimal page replacement algorithm
					work item ordering loop or no */

int    globalNextWorker = 0;
bool   stopAnnealing    = false;
bool   firstAssignment  = true;

/**
 ** END DISTRIBUTED VARS
 **************************************************
 **************************************************/



int doLinkingRecurseVisits = 0, buildTracksAddToResultsVisits = 0, compatibleEndpointsFound = 0;
#define RACE_TO_MAX_COMPATIBLE false
#define MAX_COMPATIBLE_TO_FIND 250000

int rejectedOnVelocity, rejectedOnPosition, wereCompatible;



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




/***************************************************************************
 * HELPER CLASSES
 * 
 * These should be identical to those in sequential linkTracklets.
 ***************************************************************************/


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
  TreeNodeAndTime(){
  }
  TreeNodeAndTime(lsst::mops::TrackletTreeNode * tree, ImageTime i) {
    myTree = tree;
    myTime = i;
  }
  lsst::mops::TrackletTreeNode * myTree;
  ImageTime myTime;
};




/**********************************************************************
 * 
 *  Helper structs from matt cleveland, used for distributed processing
 *
 ***********************************************************************/


/***************************************************************************
 * This is the data that is stored in the node cache.  A node cache entry
 * contains the unique identifiers for the TreeNodeAndTime object (treeId,
 * nodeId) and the TreeNodeAndTime object itself.
 ***************************************************************************/
typedef struct {
  unsigned int treeId, nodeId;
  TreeNodeAndTime node;
  bool inUse;
  time_t timeStamp;
} cachedTreeNode;

/*******************************************************************************
 * This struct holds the arguments sent to the killProcs function.  The sentinel
 * value indicates to the workers what action they should take after processing
 * their current work item set.
 *******************************************************************************/
typedef struct {
  int numProcs;
  int sentinel;
} killSentinel;

/******************************************************************************
 * This struct is used for keeping track of a work item and meta data about it, 
 * which is used when the work items are distributed.
 ******************************************************************************/
typedef struct {
  std::vector<std::pair<unsigned int, unsigned int> > endpoints;
  int numCompatibleEndpoints;
  int numNodes;
  unsigned long long timeUnits;
} workItemFile;

/*******************************************************************************
 * An index specifying a tree, by its ID, and a node in that tree, by its ID
 *******************************************************************************/
typedef struct {
  unsigned int treeId;
  unsigned int nodeId;
} treeIdNodeIdPair;






/***************************************************************************
 *
 * Prototypes for standard linkTracklets functions - jmyers
 * 
 ***************************************************************************/

void setTrackletVelocities(
    const std::vector<MopsDetection> &allDetections,
    std::vector<Tracklet> &queryTracklets);

void makeTrackletTimeToTreeMap(
    const std::vector<MopsDetection> &allDetections,
    std::vector<Tracklet> &queryTracklets,
    std::map<ImageTime, TrackletTree > &newMap,
    const linkTrackletsConfig &myConf);

unsigned int buildTracksAddToResults(
    const std::vector<MopsDetection> &allDetections,
    const std::vector<Tracklet> &allTracklets,
    linkTrackletsConfig searchConfig,
    TreeNodeAndTime &firstEndpoint,
    TreeNodeAndTime &secondEndpoint,
    const std::vector<treeIdNodeIdPair> &supportNodes,
    TrackSet & results,
    int &cacheSize, cachedTreeNode *nodeCache,
    std::map<ImageTime, TrackletTree > &trackletTimeToTreeMap,
    const std::vector<std::vector<treeIdNodeIdPair> > &finalEndpointOrder,
    unsigned long long &pageFaults, unsigned int &myNumCompatible,
    unsigned int setNum);

void recenterDetections(std::vector<MopsDetection> &allDetections, 
                        const linkTrackletsConfig &searchConfig);


bool endpointTrackletsAreCompatible(
    const std::vector<MopsDetection> & allDetections, 
    const Track &newTrack,
    const linkTrackletsConfig &searchConfig);




/***************************************************************************
 ***************************************************************************
 **  distributed linkTracklets helper functions
 **    -MGC
 ***************************************************************************
 ***************************************************************************/



void waitForTask(int rank,
		 const std::vector<lsst::mops::MopsDetection> &allDetections, //from MAIN
		 std::vector<lsst::mops::Tracklet> &allTracklets, //from MAIN
		 linkTrackletsConfig searchConfig  /*from MAIN*/);
	



/***************************************************************************
 * DLT function prototypes
 ***************************************************************************/
void 
distributeCurrentWorkload(std::vector<std::vector<workItemFile> > &assignment);

void 
*killProcs(killSentinel, const std::vector<MopsDetection> & allDets);

int
leastRecentlyUsed(cachedTreeNode *nodeCache, int cacheSize);

int 
findFarthestNodeIndex(cachedTreeNode *nodeCache, int cacheSize,
		      const std::vector<std::vector<treeIdNodeIdPair> > &finalEndpointOrder,
		      unsigned int setStart, unsigned int workItemStart);

lsst::mops::TrackletTree *
findTreeById(const unsigned int id, 
	     std::map<ImageTime, lsst::mops::TrackletTree > &treeMap);

ImageTime 
findImageTimeForTree(const unsigned int id, 
                     const std::map<ImageTime, lsst::mops::TrackletTree > &treeMap);

lsst::mops::TrackletTreeNode *
getNodeByIDAndTime(lsst::mops::TrackletTree *myTree, unsigned int nodeId);

int 
loadNodeFromCache(unsigned int treeId, unsigned int nodeId, int &cacheSize, cachedTreeNode *nodeCache,
		  std::map<ImageTime, lsst::mops::KDTree <unsigned int> > &trackletTimeToTreeMap,
		  const std::vector<std::vector<treeIdNodeIdPair> > &finalEndpointOrder,
		  unsigned long long &pageFaults);


/***************************************************
 * From http://www.parashift.com/c++-faq-lite/misc-technical-issues.html#faq-39.1
 ***************************************************/
std::string stringify(size_t x)
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

/***************************************************************************
 * Equality operator for treeIdNodeIdPair struct
 ***************************************************************************/
bool sameTreeNode(const treeIdNodeIdPair &one, const treeIdNodeIdPair &two){
  if(one.nodeId == two.nodeId && one.treeId == two.treeId){
    return true;
  }
  return false;
}

/***************************************************************************
 * Create and return time string for logging output
 ***************************************************************************/
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





void printCacheAndWorkItems(cachedTreeNode *nodeCache,
			    const std::vector<std::vector<treeIdNodeIdPair> > &finalEndpointOrder,
			    int e1, int e2, int index)
{
  
  std::ofstream outFile;
  std::string writeName = "cache_analysis_" +  stringify(rand() % 10000000);
  outFile.open(writeName.c_str());

  
  outFile << "Looking for treeId/nodeId: " << e1 << "/" << e2 << std::endl;

  if(index != -1){
    outFile << "Found pair at index: " << index << std::endl;
  }
  outFile << "------------------------------" << std::endl << std::endl;
  

  outFile << "------------------------------" << std::endl;
  outFile << "|   node cache contents      |" << std::endl;
  outFile << "------------------------------" << std::endl;

  for(unsigned int i=0; i < globalCacheSize; ++i){
    outFile << "NodeCache [" << i << "]: " << nodeCache[i].treeId << "/" << nodeCache[i].nodeId << std::endl;
  }
  outFile << "|                            |" << std::endl;
  outFile << "------------------------------" << std::endl;

  outFile << std::endl << "------------------------------" << std::endl;
  outFile <<              "|  work items                |" << std::endl;
  outFile << "------------------------------" << std::endl;
  for(unsigned int i=0; i<finalEndpointOrder.size(); ++i){
    outFile << "**Work item " << i << std::endl;
    for(unsigned int j=0; j<finalEndpointOrder.at(i).size(); ++j){
      outFile << finalEndpointOrder.at(i).at(j).treeId << "/"
		<< finalEndpointOrder.at(i).at(j).nodeId << std::endl;
    }
  }
  outFile.close();
}



/***************************************************************************
 * Calculate size of treeNode in bytes.  Read that number of bytes from a 
 * random file and then return.  This is used to demonstrate the cost of 
 * disk accesses when cache misses occur.
 ***************************************************************************/
#define FILE_SIZE 10000
unsigned long 
readNodeFromDisk(TrackletTreeNode *treeNode)
{

  unsigned long numPointsAndValues = 0;
  for(unsigned int i=0; i < treeNode->getMyData()->size(); ++i){
    numPointsAndValues += treeNode->getMyData()->at(i).getPoint().size();
  }

  unsigned long numBytes = ((sizeof(unsigned int) * 4) +
			    (sizeof(double) * treeNode->getUBounds()->size()) +
			    (sizeof(double) * treeNode->getLBounds()->size()) +
			    (sizeof(unsigned int)) +
			    (sizeof(double) * numPointsAndValues));
  
  
  int offset = rand() % (FILE_SIZE - numBytes);

  FILE *readFile;
  std::string fileName = "diskSimulation/disk_read_file_" + stringify(rand() % FILE_SIZE);
  readFile = fopen(fileName.c_str(), "r");
  fseek(readFile, offset, SEEK_SET);
  
  unsigned long bytesRead = 0;
  int c;

  if(readFile == NULL){
    std::cerr << timestamp() << "Unable to open file: " << fileName << std::endl;
    return -1;
  }
  else{
    do{
      c = fgetc(readFile);
      ++bytesRead;
    }while(c != EOF && bytesRead < numBytes);
    fclose(readFile);
  }
  return numBytes;
}


/***************************************************************************
 * Check the cache for the requested node, if it's there, return it, if not
 * load it into the cache using the optimal page replacement algorithm.
 *   -MGC
 ***************************************************************************/
int loadNodeFromCache(unsigned int treeId, unsigned int nodeId, int &cacheSize, cachedTreeNode *nodeCache,
		      std::map<ImageTime, TrackletTree > &trackletTimeToTreeMap,
		      //const std::vector<std::vector<int> > &finalEndpointOrder,
		      const std::vector<std::vector<treeIdNodeIdPair> > &finalEndpointOrder,
		      unsigned long long &pageFaults, unsigned int &currPf,
		      unsigned int setNum, unsigned int nodeNum)
{
  bool fullCache = true;

  //if cache has free slots, load the node in automatically
  if(globalCacheSize < CACHE_SIZE){
    fullCache = false;
  }

  //if this TAT is already in the cache, load it and skip the OPR section
  for(unsigned int i=0; i < globalCacheSize; ++i){
    if(nodeCache[i].nodeId == nodeId && nodeCache[i].treeId == treeId){
      ++globalCacheHits;
      nodeCache[i].timeStamp = std::clock();
      return i;
    }
  }

  ++globalCacheMisses;
  ++currPf;
  ++pageFaults;

  //The TAT was not found in the cach, create new TAT.  Place this new TAT in the
  //cache, evicting the TAT that will be used the farthest in the future
  int mostDistantIndex;

  //if the cache is full, use OPR to find the index to load this node into
  if(fullCache){
    
    //this is the index in the cache of the node that will be used 
    //farthest in the future, the basis for optimal page replacemtn
    if(DO_OPRA){
      mostDistantIndex = findFarthestNodeIndex(nodeCache, cacheSize, finalEndpointOrder,
					       setNum, nodeNum);
    }
    //use the Least Recently Used cache replacment algorithm
    else if(DO_LRU){
      mostDistantIndex = leastRecentlyUsed(nodeCache, cacheSize);
      if(mostDistantIndex == -1){
	std::cerr << timestamp() << "got -1 from LRU" << std::endl;
      }
    }
    //hoose random node to evict    
    else{
      mostDistantIndex = rand() % CACHE_SIZE;

      //make sure aren't using a used node
      while(nodeCache[mostDistantIndex].inUse){
	mostDistantIndex = rand() % CACHE_SIZE;
      }
    }
  }
  //if cache is not full, load the node into the next available cache slot
  else{
    mostDistantIndex = globalCacheSize;
    ++globalCacheSize;
  }
  
  //Get the tree and ImageTime for this image (tree) ID
  TrackletTree *myTree = findTreeById(treeId, trackletTimeToTreeMap);
  ImageTime it = findImageTimeForTree(treeId, trackletTimeToTreeMap);
  
  //find this node in its tree
  TrackletTreeNode *treeNode = NULL;
  if( myTree != NULL ){
    treeNode = getNodeByIDAndTime(myTree, nodeId);
  }
  
  if( treeNode != NULL ){
    /**********************************************************************
     * Simulate reading from disk to enable accurate timing data for
     * performance analysis
     **********************************************************************/
    if(SIMULATE_DISK_ACCESS)
      readNodeFromDisk(treeNode);

    nodeCache[mostDistantIndex].nodeId = nodeId;
    nodeCache[mostDistantIndex].treeId = treeId;

    nodeCache[mostDistantIndex].node.myTree = treeNode;
    nodeCache[mostDistantIndex].node.myTime = it;
    
    nodeCache[mostDistantIndex].timeStamp = std::clock();

    return mostDistantIndex;
  }
  else{
    std::cerr << timestamp() << "TAT not found, node cache unchanged." << std::endl;
  }

  return (-1);
}



/**
 ** END Helper functions
 ***************************************************************************
 ***************************************************************************/




/***************************************************************************
 ***************************************************************************
 **  DISTRIBUTED FUNCTION DEFINITIONS
 **    -MGC
 ***************************************************************************
 ***************************************************************************/

/***************************************************************************
 * Given an ImageTime id, find the corresponding tree and return it.
 ***************************************************************************/
TrackletTree *
findTreeById(const unsigned int id, 
	     std::map<ImageTime, TrackletTree > &treeMap)
{
  TrackletTree *retVal = NULL;

  std::map<ImageTime, TrackletTree >::iterator mapIter;
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
			       const std::map<ImageTime, TrackletTree > &treeMap)
{

  std::map<ImageTime, TrackletTree >::const_iterator mapIter;
  //MGC 1/30/11
  ImageTime it;
  for(mapIter = treeMap.begin(); mapIter != treeMap.end(); mapIter++){
    it = (*mapIter).first;
    if(it.getImageId() == id){
      return it;
    }
  }

  ImageTime retVal;
  return retVal;
}



TrackletTreeNode *
getNodeByIDAndTime(TrackletTree *myTree, unsigned int nodeId)
{
  TrackletTreeNode *retVal;

  //search tree for nodeId
  std::vector<TrackletTreeNode *> nodeList;
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



/**********************************************************************
 * Matt using so each processor can write out its own results
 *   -MGC
 **********************************************************************/
void writeMyResults(std::string outFileName, 
		    const std::vector<MopsDetection> &allDets,
		    const std::vector<Tracklet> &allTracklets,
		    TrackSet &toRet, unsigned long long nextNo)
{

  // create the output directory for this worker, if it doesn't exist
  if (nextNo == 0) {
    std::string newPath = "results/" + outFileName;
    std::cerr << "Making directory " << outFileName;
    int rc = mkdir(newPath.c_str(), 0777);
    if (rc) {
      std::cerr << "Failed to create " + newPath + " directory with " << rc << std::endl;
      exit(rc);
    }
  }

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
void collectTrackResults(long numProcs, const std::vector<MopsDetection> & allDets)
{
  MPI_Status status;

  int _from;
  for(_from = 1; _from < numProcs; ++_from){

    //determine number of tracks to expect
    unsigned int numTracks;
    MPI_Recv(&numTracks, 1, MPI_INT, _from, MPI_COLLECT_TAG, MPI_COMM_WORLD, &status); 

    if (numTracks > 0) {
      for(unsigned int i=0; i < numTracks; ++i){
	
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
	Track thisTrack;
	
	thisTrack.componentTrackletIndices  = tSet;
	
	std::set<unsigned int>::const_iterator detIter;
	for (detIter  = dSet.begin();
	     detIter != dSet.end();
	     detIter++){
	  thisTrack.addDetection(*detIter, allDets);
	}
	
	toRet->insert(thisTrack);
      }
    }
  }
  std::cerr << timestamp() << "Master exiting collectTracks with toRet size " << toRet->size() << std::endl;
  //stopAnnealing = true; //matt added 3/3/11  4/1/11
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
  unsigned int ti = -1;
  unsigned int ni = -1;
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
 * Determine the cost of this set work item the exeuction model.
 ***************************************************************************/
unsigned long long executionModelCost(workItemFile &myItem)
{
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
void generateNeighbor(std::vector<std::vector<workItemFile> > &assignment, 
		      std::vector<std::vector<workItemFile> > &retVal)
{
  std::vector<std::vector<workItemFile> > mixedSet;

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

  //divvy out the work items
  unsigned long long timePerWorker  = (totalTimeUnits / numWorkers);
  unsigned int numAssignedWorkItems = 0;
  unsigned long long assignedTime   = 0;
  
  std::vector<workItemFile> temp;
  bool full = false;

  for(unsigned int i=0; i < mixedSet.size() && !full; ++i){
    for(unsigned int j=0; j < mixedSet.at(i).size() && !full; ++j){

      temp.push_back(mixedSet.at(i).at(j));
      assignedTime += executionModelCost(mixedSet.at(i).at(j));
      ++numAssignedWorkItems;

      if(assignedTime >= timePerWorker){
	retVal.push_back(temp);
	temp.clear();
	assignedTime = 0;

	if(retVal.size() == mixedSet.size()){
	  full = true;
	}
      }
    }
  }

  //assign any leftover work items
  unsigned int numRemaining = totalWorkItems - numAssignedWorkItems;
  if(numRemaining > 0){

    unsigned int nextWorker = 0;
    
    unsigned int n = mixedSet.size() - 1;
    int remainder = numRemaining - mixedSet.at(n).size();

    while(remainder > 0){
      --n;
      remainder = remainder - mixedSet.at(n).size();
    }
    
    //this is the index into the nth work item set we should start with
    if(n != (mixedSet.size() - 1)){
      remainder = (-remainder);
      //++n;
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
unsigned long long calculateSetAssignmentCost(std::vector<std::vector<workItemFile> > &assignment)
{
  //L is the constant value representing the cost to load a node from memory
  unsigned int L = 1000;

  //this is the max cost from this work item set assignment
  unsigned long long maxCost = 0;

  for(unsigned int i=0; i < assignment.size(); ++i){
    
    unsigned long long cost = 0;

    //look through all the nodes of the adjacent work item and find number of shared nodes
    for(unsigned int j=0; j < assignment.at(i).size(); ++j){
      cost += executionModelCost(assignment.at(i).at(j));
      
      workItemFile currentWorkItem = assignment.at(i).at(j);
      unsigned int cacheHits = 0;
      unsigned int numNodes = currentWorkItem.numNodes;;
      
      /* 
       *memory (IO) model cost (this only goes here because the IO
       * cost varies at the WORK ITEM SET ASSIGNMENT level, not the work
       * item level 
       */

      //the first work item will face an empty cache, it will need to load all of its
      //nodes into memory
      if(j == 0){
	cost+=(numNodes * L);
      }
      else{
	workItemFile previousWorkItem = assignment.at(i).at((j-1));

	//only the first CACHE_SIZE nodes of this work item will be match candidates
	for(unsigned int k=0; ((k <  currentWorkItem.endpoints.size()) && (k < CACHE_SIZE)); ++k){
	  //for(unsigned int k=0; k <  currentWorkItem.endpoints.size(); ++k){

	  //only the last CACHE_SIZE nodes from the previous work item will
	  //remain in the cache, only iterate over those
	  unsigned int startIndex = previousWorkItem.endpoints.size() - CACHE_SIZE;
	  if(startIndex < 0){
	    startIndex = 0;
	  }

	  for(unsigned int l=startIndex; l < previousWorkItem.endpoints.size(); ++l){
	    if((currentWorkItem.endpoints.at(k).first == previousWorkItem.endpoints.at(l).first) &&
	       (currentWorkItem.endpoints.at(k).second == previousWorkItem.endpoints.at(l).second)){
	      ++cacheHits;
	      break;
	    }
	  }
	}
     
	cost+=(L * (numNodes - cacheHits));

      } /*end else*/
    } /* end j */

    //keep track of the greatest cost encountered
    if(cost > maxCost){
      maxCost = cost;
    }
    
  } /* end i */

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
void determineWorkItemSets(std::vector<std::vector<workItemFile> > &assignment)
{
  unsigned int numWorkers = assignment.size();

  //Vector (per worker) of vectors of work items.  This is a work item
  //set assignment
  std::vector<std::vector<workItemFile> > currentAss;
  for(unsigned int i=0; i < assignment.size(); ++i){
    currentAss.push_back(assignment.at(i));
  }

  //number of iterations to consider different work item set 
  //assignment permutations
  int maxIterations = 10;
  int iteration     = 0;
  unsigned int T    = 10;

  //record-keeping variables for the simulated annealing algorithm
  unsigned long long bestCost = ULONG_MAX;

  //ensure at least one round of "intelligent" work item 
  //distribution is considered
  bool singleRun = false;
  if(stopAnnealing){
    singleRun = true;
  }

  while(!stopAnnealing || singleRun){

    //annealing sections
    if(DO_ANNEALING){

      //this is the work item set assignment created as a result of 
      //shuffling `current'
      std::vector<std::vector<workItemFile> > trialAss;
   
      generateNeighbor(currentAss, trialAss);
      
      unsigned long long trialCost = ULONG_MAX;
      
      if(trialAss.size() == numWorkers){             //TODO REMOVE CONDITION
	trialCost = calculateSetAssignmentCost(trialAss);
      }

      //this is the new best set assignment
      if( trialCost < bestCost ){
      /* 	std::cerr << timestamp() << "updating best work item set assignment in SA. Trial cost is " 
                << trialCost << " bestCost is " << bestCost << std::endl; */
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
      if((updateCurrent(trialCost, bestCost, T)) && (trialAss.size() == numWorkers)){ //TODO REMOVE second condition
	
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
      break;
    }

    //"temperature" for P(e, e', T) function, never < 0
    if(T > 0){
      --T;
    }

    //iterations
    ++iteration;

    //sleep if not performing annealing to prevent spinning
    if(!DO_ANNEALING){
      sleep(5);
    }
  }
  
  //reset the annealing sentinel
  stopAnnealing = false;
}

struct annealArg {
  int numProcs;
};


/***************************************************************************
 * Wait to hear from all workers that they have completed processing their
 * work items.  When have heard from all, tell simulated annealing algorithm
 * to stop.
 ****************************************************************************/
void *annealingSentinel(void *arg)
{
  struct annealArg * sArg = (annealArg *)arg;
  int np = sArg->numProcs;

  std::cerr << timestamp() << "Master triggering annealing barrier, waiting on numworkers " << np << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);

  //detach this thread so its resources are freed upon completion
  pthread_detach(pthread_self());

  //set the global annealing sentinel
  stopAnnealing = true;

  std::cerr << timestamp() << "Annealing thread is dying" << std::endl;
  
  return NULL;
}



/************************************************************
 * Distribute work items to worker nodes according to the 
 * distribution defined in "workItemSets".
 ************************************************************/
void distributeWorkload(std::vector<std::vector<workItemFile> > &workItemSets)
{
  unsigned int numProcs = workItemSets.size();

  long int workUnitsList[numProcs];
  for(unsigned int i=0; i<numProcs; ++i){
    workUnitsList[i] = 0;
  }

  std::cerr << "Master distributing workload..." << std::endl;

  for(unsigned int i=0; i < workItemSets.size(); ++i){

    //corner case where a worker wasn't assigned any work items, inform him that he 
    //should expect zero items
    if(workItemSets.at(i).size() == 0){
      int numEndpoints = 0;
      MPI_Send(&numEndpoints, 1, MPI_INT, (i+1), /*MPI_NODE_COMM_TAG*/i, MPI_COMM_WORLD);
    }
    else{
      for(unsigned int j=0; j < workItemSets.at(i).size(); ++j){

	workItemFile workItem = workItemSets.at(i).at(j);
	
	//send number of endpoints to recipient
	unsigned int numEndpoints = workItem.endpoints.size();

	//send to (i+1), i starts at 0 but proc 0 doesn't process work items
	MPI_Send(&numEndpoints, 1, MPI_INT, (i+1), /*MPI_NODE_COMM_TAG*/j, MPI_COMM_WORLD);
	
	//keep track of predicted time required at each processor
	workUnitsList[i] += workItem.timeUnits;
	
	//aggregate data to send from each endpoint
	int bufferSize = numEndpoints * 2;
	int sendBuffer[bufferSize];
	int index = 0;
	
	for(unsigned int k=0; k < numEndpoints; ++k){
	  
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
	MPI_Send(&sendBuffer, bufferSize, MPI_INT, (i+1), /*MPI_NODE_COMM_TAG*/j, MPI_COMM_WORLD);
      }
    }
  } /* end file transmission */


  //just for debug
  for(unsigned int i=1; i <= numProcs; ++i){
    std::string myTime("Processor ");
    /* myTime = timestamp() + myTime + stringify(i) + " got " + stringify(workUnitsList[(i-1)]) 
       + " time units and " + stringify(workItemSets.at((i-1)).size()) + " total work items.\n"; 
       std::cerr << myTime; */
  }

  //this thread waits to hear back from the workers that they have completed
  //processing on the work items they were assigned
  struct annealArg arg;
  arg.numProcs = numProcs;
  int rc = pthread_create(&thread, NULL, annealingSentinel, (void*)&arg);
  std::cerr << timestamp() << "Master creating thread with numProcs: " << arg.numProcs << std::endl;
  if(rc){
    std::cerr << "pthread_create return code is" << rc << ", aborting." << std::endl;
    exit(rc);
  }
}



/**********************************************************************
 * Tell all the slave nodes they can stop working.
 **********************************************************************/
void *killProcs(killSentinel ta, const std::vector<MopsDetection> & allDets)
{
  int _sentinel = ta.sentinel;
  
  for(int i=1; i < ta.numProcs; ++i){
    MPI_Send(&_sentinel, 1, MPI_INT, i, i, MPI_COMM_WORLD);
  }

  //collect results and add to return vector if workers not writing to disk
  if(!WRITE_RESULTS_TO_DISK){
    collectTrackResults(ta.numProcs, allDets);
  }

  return NULL;
}



/***************************************************************************
 * Find timestamp with smallest value
 ***************************************************************************/
int leastRecentlyUsed(cachedTreeNode *nodeCache, int cacheSize)
{
  int LRU = -1;
  
  double oldest = std::numeric_limits<double>::max();
  
  std::cerr.setf(std::ios::fixed);

  for(unsigned int i=0; i < globalCacheSize; ++i){

    if(nodeCache[i].timeStamp < oldest && !nodeCache[i].inUse){
      oldest = nodeCache[i].timeStamp;
      LRU = i;
    }
  }

  return LRU;
}



/***************************************************************************
 * Look through the cache and see which cached node is used farthest in
 * the future.  If a cached node is found that is never used again, there
 * is no need to check any more nodes.
 ***************************************************************************/
int findFarthestNodeIndex(cachedTreeNode *nodeCache, int cacheSize,
			  const std::vector<std::vector<treeIdNodeIdPair> > &finalEndpointOrder,
			  unsigned int setStart, unsigned int nodeStart)
{
  int farthest = -1;
  unsigned int farthestOffset = 0;

  //calculate total number of nodes in this work item set assignment
  unsigned int totalNumNodes = 0;
  for(unsigned int i=0; i < finalEndpointOrder.size(); ++i){
    totalNumNodes += finalEndpointOrder.at(i).size();
  }

  for(unsigned int i=0; i < globalCacheSize; ++i){

    if(!nodeCache[i].inUse){
  
      unsigned int currentOffset = totalNumNodes;
      bool everFound = false;
      
      //iterate over the work items to see where the ith entry in the node cache is last used
      std::vector<std::vector<treeIdNodeIdPair> >::const_reverse_iterator setIterator;
      std::vector<treeIdNodeIdPair>::const_reverse_iterator workItemIterator;

      unsigned int setCount = finalEndpointOrder.size();
      
      for(setIterator = finalEndpointOrder.rbegin();
	  setIterator != finalEndpointOrder.rend() && setCount >= setStart;
	  ++setIterator){

	unsigned int nodeCount = (*setIterator).size();

	for(workItemIterator = (*setIterator).rbegin();
	    workItemIterator != (*setIterator).rend() && 
	      (setCount >= setStart && nodeCount > nodeStart);
	    ++workItemIterator){

	  //check if this is the same node
	  if((*workItemIterator).treeId == nodeCache[i].treeId &&
	     (*workItemIterator).nodeId == nodeCache[i].nodeId){

	    if(currentOffset > farthestOffset){
	      farthestOffset = currentOffset;
	      farthest = i;
	    }
	    everFound = true;
	  }
	  
	  --currentOffset;
	  --nodeCount;
	} // end work item node check
	
	--setCount;
      } // end work item set check

      //node not found, therefore never used again
      if(!everFound){
	return i;
      }
    } // in use check 
  } // end nodeCache iterations (i) 

  return farthest;
}






/************************************************************
 * Slave processors wait for signals from master.
 * They accept the data and process it in doLinkingRecurse2.
 ************************************************************/
void waitForTask(int rank,
		 std::vector<MopsDetection> &allDetections,
		 std::vector<Tracklet> &allTracklets, 
		 linkTrackletsConfig searchConfig)
{
  //this is the set of tracks this function returns
  TrackSet myTracks;

  //this is THE node cache, declare and initialize it
  cachedTreeNode nodeCache[CACHE_SIZE];
  for(unsigned int i=0; i < CACHE_SIZE; ++i){
    nodeCache[i].nodeId = -1;
    nodeCache[i].treeId = -1;
    nodeCache[i].inUse = false;
    nodeCache[i].timeStamp = std::clock();
  }

  //number of elements in the cache
  int cacheSize = 0;
  

  //need to make the time/tree map so we can identify common Node IDs and 
  //therefore reduce master/slave communication to only those IDs that 
  //require linking

  // jmyers: we now need to do recentering just like in the linkTracklets() function
  recenterDetections(allDetections, searchConfig);
  setTrackletVelocities(allDetections, allTracklets);

  std::map<ImageTime, TrackletTree > trackletTimeToTreeMap;
  makeTrackletTimeToTreeMap(allDetections, allTracklets, 
                            trackletTimeToTreeMap, searchConfig);
  std::cerr << timestamp() << "Worker " << rank << " made trees" << std::endl;

  //timing stuff
  int numWorkItems = 0;
  double myBFTime = 0;

  std::vector<std::vector<treeIdNodeIdPair> > allEndpoints;
  unsigned long long myTotalTracks = 0;

  //this is for debugging
  std::ofstream outFile;
  std::string writeName = "rank/" +  stringify(rank);
  outFile.open(writeName.c_str());
  std::cerr << "Worker " << rank << " entering infinite loop..." << std::endl;

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
    MPI_Recv(&numEndpoints, 1, MPI_INT, 0/*MPI_ANY_SOURCE*/,
	     /* MPI_NODE_COMM_TAG*/MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    /* std::cerr << timestamp() << "Worker " << rank << " got " << numEndpoints << " endpoints" << std::endl; */

    /***************************************************************************
     * numEndpoints == -1 is master node's signal to stop working
     * numEndpoints == -2 indicates a work item set assignment is complete
     * numEndpoints == 0 is an empty work item set assignment
     ***************************************************************************/
    if (numEndpoints == 0 || numEndpoints == -1 || numEndpoints == -2) {
        /* std::cerr << timestamp() << "Worker " << rank << " got " << numEndpoints << " endpoints" << std::endl; */
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
      if (allEndpoints.empty() || numEndpoints == 0) {
	
	//this empty endpoint set happens to be the last work item set assignment
	//corner case discovered in testing
	if(numEndpoints == -1){// || stopAnnealing){  //TODO WHAT IS THIS SECOND CONDITION HERE FOR?
	  outFile.close();

	  //matt added this on 3/6/11
	  if(!WRITE_RESULTS_TO_DISK){
	    int writing = 0;
	    std::cerr << timestamp() << "Worker " << rank << " sending termination signal to master" << std::endl;
	    MPI_Send(&writing, 1, MPI_INT, 0, MPI_COLLECT_TAG, MPI_COMM_WORLD);
	  }
	  std::cerr << timestamp() << "Worker " << rank << " triggering early barrier" << std::endl;

	  MPI_Barrier(MPI_COMM_WORLD); //matt 3/6/11
	  std::cerr << timestamp() << "Worker " << rank << " has passed the early trigger barrier" << std::endl;
	  break;
	}
	else{
	  //tell the master we are finished and ready for next assignment
	  std::cerr << timestamp() << "Worker " << rank << " waiting on barrier at line " << __LINE__ << std::endl;
	  //MPI_Barrier(MPI_COMM_WORLD); //matt 3/6/11 LOOK AT THIS 3/30/11
	  continue;
	}
      }

      std::vector<std::vector<treeIdNodeIdPair> > finalEndpointOrder;

      //use allEndpoints to iterate through the work items and 
      //sort them accordingly
      double opra_1_end = std::clock();
      if(DO_OPRA_PREP){
	double opra_1_start = std::clock();

	finalEndpointOrder.push_back(allEndpoints.at(0));
	allEndpoints.erase(allEndpoints.begin());
	
	while(allEndpoints.size() > 0){
	  
	  //tally the number of common endpoints of each work item
	  int matches[allEndpoints.size()];

	  //initialize tally
	  for(unsigned int m=0; m < allEndpoints.size(); ++m){
	    matches[m] = 0;
	  }
	  
	  unsigned int latestWorkItemIndex = (finalEndpointOrder.size() - 1);
	  unsigned int latestWorkItemSize = finalEndpointOrder.at(latestWorkItemIndex).size();
	  
	  int bestMatch = 0;
	  int matchIndex = 0;

	  //only the first CACHE_SIZE nodes of this work item will be match candidates
	  unsigned int startIndex = latestWorkItemSize - CACHE_SIZE;
	  if(startIndex < 0){
	    startIndex = 0;
	  }
	  
	  for(unsigned int i=startIndex; ((i < latestWorkItemSize) && (i < CACHE_SIZE)); i+=2){
	    
	    //compare this work item to all the remaining work items
	    for(unsigned int k=0; k < allEndpoints.size(); ++k){

	      //only the last CACHE_SIZE nodes from the previous work item will
	      //remain in the cache, only iterate over those
	      for(unsigned int l=0; ((l < allEndpoints.at(k).size()) && (l < CACHE_SIZE)); ++l){
		if( sameTreeNode(finalEndpointOrder.at(latestWorkItemIndex).at(i), allEndpoints.at(k).at(l)) ){
		  int temp = matches[k];
		  ++temp;
		  matches[k] = temp;
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

	opra_1_end = getTimeElapsed(opra_1_start);
      }
      //do not do the OPR work item massaging, process the work items in the order in 
      //which they were received
      else{
	for(unsigned int i=0; i < allEndpoints.size(); ++i){
	  finalEndpointOrder.push_back(allEndpoints.at(i));
	}
	for(unsigned int i=0; i < allEndpoints.size(); ++i){
	  allEndpoints.at(i).clear();
	}
	allEndpoints.clear();
      }
      

      //iterate through finalEndpointOrder vector and perform "track extraction"
      //on its work item
      //NOTE: finalEndpointOrder is a vector of work items
      double totalOpra2Time = 0;
      unsigned long long pageFaults = 0;
      for(unsigned int i=0; i < finalEndpointOrder.size(); ++i){
	
	double opra_2_start = std::clock();

	//set up endpoint values for brute force computation
	int e1 = -1, e2 = -1;
	std::vector<treeIdNodeIdPair> supportPoints;
	unsigned int treeId1 = 0;
	unsigned int nodeId1 = 0; 
	unsigned int treeId2 = 0;
	unsigned int nodeId2 = 0;

	for(unsigned int j=0; j < finalEndpointOrder.at(i).size(); ++j){

	  //first iteration receives first endpoint
	  unsigned int currPf = 0;
	  if(j == 0){
	    treeId1 = finalEndpointOrder.at(i).at(j).treeId;
	    nodeId1 = finalEndpointOrder.at(i).at(j).nodeId;
	    e1 = loadNodeFromCache(treeId1, nodeId1, cacheSize, 
				   &nodeCache[0], trackletTimeToTreeMap,
				   finalEndpointOrder, pageFaults, currPf, i, 1);

	    //invalid node cache lookup, abort this work item
	    if(e1 == -1){
	      std::cerr << timestamp() << "Firstendpoint is NULL" << std::endl;
	      break;
	    }
	    else{
	      nodeCache[e1].inUse = true;
	    }
	  }
	  //second iteration receives second endpoint
	  else if(j == 1){
	    treeId2 = finalEndpointOrder.at(i).at(j).treeId;
	    nodeId2 = finalEndpointOrder.at(i).at(j).nodeId;
	    e2 = loadNodeFromCache(treeId2, nodeId2, cacheSize, 
				   &nodeCache[0], trackletTimeToTreeMap,
				   finalEndpointOrder, pageFaults, currPf, i, 2);

	    //invalid node cache lookup, abort this work item
	    if(e2 == -1){
	      std::cerr << timestamp() << "Second endpoint is NULL" << std::endl;
	      break;
	    }
	    else{
	      nodeCache[e2].inUse = true;
	    }
	  }
	  //all other iterations are support points
	  else{
	    supportPoints.push_back(finalEndpointOrder.at(i).at(j));
	  }
	}

	double opra_2_end = getTimeElapsed(opra_2_start);
	totalOpra2Time += opra_2_end;
	double workItemStart = 0, wTime = 0;

	//sanity check
	if(e1 == -1 || e2 == -2){
	  std::cerr << timestamp() << "Endpoint cache miss occurred, e1 == " << e1 << ", e2 == " << e2 << std::endl;
	}
	//further sanity check
	else if(e1 == e2){
	  std::cerr << timestamp() << "Rank " << rank << " has an invalid 1st/2nd endpoint pair.  They are the same node.  Aborting." << std::endl;
	  std::cerr << timestamp() << "Rank " << rank << " first from cache: " << nodeCache[e1].treeId << " / " << nodeCache[e1].nodeId 
		    << " at index " << e1 << " looking for " << treeId1 << " / " << nodeId1 << std::endl;
	  std::cerr << timestamp() << "Rank " << rank << " second from cache: " << nodeCache[e2].treeId << " / " << nodeCache[e2].nodeId 
		    << " at index " << e2 << " looking for " << treeId2 << " / " << nodeId2 << std::endl;
	}
	//build the track
	else if(e1 != -1 && e2 != -1){

	  workItemStart = std::clock();
	  unsigned int myNumCompatible = 0;
	  unsigned int pf = buildTracksAddToResults(allDetections, allTracklets, searchConfig,
						    nodeCache[e1].node, nodeCache[e2].node,
						    supportPoints, 
						    myTracks, cacheSize, &nodeCache[0],
						    trackletTimeToTreeMap,
						    finalEndpointOrder, pageFaults, myNumCompatible, i);

	  wTime = getTimeElapsed(workItemStart);
	  outFile << timestamp() << "Worker " << rank << " built track in " << wTime << " seconds, had " 
		  << pf << " page faults, had " << supportPoints.size() << " support points and "
		  << myNumCompatible << " compatible tracklets." << std::endl;
	}
	else{
	  std::cerr << "Cache corruption has occurred while loading one of the endpoints." << std::endl;
	}
	
	//reset node cache use flags
	nodeCache[e1].inUse = false;
	nodeCache[e2].inUse = false;

	//TIMING STUFF
	myBFTime += wTime;
      }
      /*
       * END OPTIMAL PAGE REPLACEMENT
       **************************************************/


      //output processing time
      double diff = getTimeElapsed(start);
      std::string thisOut("Rank ");
      outFile << timestamp() << thisOut << rank << " got " << numWorkItems << " work items" 
	      << ", spent " << diff << " seconds processing, " << opra_1_end << " seconds in OPRA 1, " 
	      << totalOpra2Time << " seconds in OPRA 2, " 
	      << myBFTime << " seconds in brute force, and found " << myTracks.size() << " tracks with " 
	      << pageFaults << " page faults and " << globalCacheHits << " cache hits." << std::endl;
      
      //reset timer
      myBFTime = 0;
      
      //not returning results to master, write them to disk
      if(WRITE_RESULTS_TO_DISK){
	writeMyResults(stringify(rank), allDetections, allTracklets, myTracks, myTotalTracks);
      }
      //return results to the master
      else{
	//return all of the Tracks found, if any
	//send any track data back to the master node
	unsigned int resultSetSize = myTracks.size();

	//inform the master how many Tracks to expect
	std::cerr << timestamp() << "worker " << rank << " responding to master with " << resultSetSize << " tracks" << std::endl;
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
	      componentTrackletIndices[count] = (int)*cIt; ++count;
	    }
	    

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
	      componentDetectionIndices[count] = (int)*dIt; ++count;
	    }

	    MPI_Send(&numDetections, 1, MPI_INT, 0, MPI_COLLECT_TAG, MPI_COMM_WORLD);
	    MPI_Send(&componentDetectionIndices, numDetections, MPI_INT, 0, MPI_COLLECT_TAG, MPI_COMM_WORLD);
	  }
	}
      }

      //tell the master we are finished and ready for next assignment
      std::cerr << timestamp() << "Worker " << rank << " triggering barrier" << std::endl;
      MPI_Barrier(MPI_COMM_WORLD);
      
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
	outFile.close();
	break;
      }
      else{
	std::cerr << timestamp() << "Unexpected number of endpoints received in termination block. Aborting." << std::endl;
	outFile.close();
	break;
      }
     
    }
    //this is a legitimate work item -- all we do here is buffer the endpoint indices of
    //all the work items we are assigned.  Once we have been assigned all of our work items,
    //we sort them and perform "track extraction" on them in that order.
    else{
      
      //receive the node index data
      int bufferSize = numEndpoints * 2;
      int indexBuffer[bufferSize];
      
      MPI_Recv(&indexBuffer, bufferSize, MPI_INT, 0, /*MPI_NODE_COMM_TAG*/MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      //update the list of endpoints
      std::vector<treeIdNodeIdPair> temp;

      unsigned int thisTreeId = 0, thisNodeId = 0, lastTreeId = 0, lastNodeId = 0;
      for(int i=0; i < bufferSize;){
	treeIdNodeIdPair p;
	
	lastTreeId = thisTreeId;
	lastNodeId = thisNodeId;

	p.treeId = indexBuffer[i];
	++i;
	
	thisTreeId = indexBuffer[i];

	p.nodeId = indexBuffer[i];
	++i;

	thisNodeId = indexBuffer[i];

	if((thisTreeId == lastTreeId) && (thisNodeId == lastNodeId)){
	  std::cerr << timestamp() << "Rank " << rank << " got sequentially identical treeId/nodeIds" << std::endl;
	}

	temp.push_back(p);

      }
      
      allEndpoints.push_back(temp);
      ++numWorkItems;
    } /* end termination check */
  } /* end while(true) loop */
}



//print the work item sets per 1/26 meeting with DKL & MGC
void printTuples(const std::vector<std::vector<workItemFile> > &assignment)
{
  std::cerr << "Printing tuples" << std::endl;

  std::string tupleFile = "tuple_comparison_";
  tupleFile = tupleFile + stringify(numGenerations);
  std::ofstream tupleOut;
  tupleOut.open(tupleFile.c_str());

  for(unsigned int i=0; i<assignment.size(); ++i){
    for(unsigned int j=0; j<assignment.at(i).size(); ++j){
      for(unsigned int k=0; k<assignment.at(i).at(j).endpoints.size(); ++k){
	tupleOut << numGenerations << ", " << i << ", " << j << ", " 
		 << assignment.at(i).at(j).endpoints.at(k).first  << ", "
		 << assignment.at(i).at(j).endpoints.at(k).second << std::endl;
      }
    }
  }

  
  tupleOut.close();

  std::cerr << "Printing tuples completed" << std::endl;
}


/**********************************************************************************
 * perform simulated annealing then send the work items to the workers
 * as calculated by that procedure
 **********************************************************************************/
void distributeCurrentWorkload(std::vector<std::vector<workItemFile> > &assignment)
{
  //only have one node distribution issued at a time, we do the simulated annealing
  //step if we are already waiting on threads to complete their processing. SA takes
  //a while, so may as well make that waiting useful.
  //determine which work items will be assigned to which workers
  double determineAssignmentStart = std::clock();
  determineWorkItemSets(assignment);
  double assTime = getTimeElapsed(determineAssignmentStart);
  std::cerr << timestamp() << "Master determined work item set assignment in " << assTime << " seconds." << std::endl;

  //master node will distribute the work items created during the recursive tree walk
  //to the worker nodes
  distributeWorkload(assignment);
}







//jmyers - added this function to increase readability of
// doLinkingRecurse.  This is the code which was written by Matt to
// add a work item, replacing our call to buildTracksAddToResults in
// doLinkingRecurse

void addWorkItem(const std::vector<MopsDetection> &allDetections,
                 const std::vector<Tracklet> &allTracklets,
                 linkTrackletsConfig searchConfig,
                 TreeNodeAndTime &firstEndpoint,
                 TreeNodeAndTime &secondEndpoint,
                 std::vector<TreeNodeAndTime> &newSupportNodes,
                 int numProcs, 
                 std::map<ImageTime, TrackletTree > &trackletTimeToTreeMap,
                 std::vector<std::vector<workItemFile> > &assignment)
{

                
    /************************************************************
     * Calculate the number of compatible endpoints in this set.
     * This allows us to predict the amount of work required by
     * this work item and distributed it accordingly.
     ************************************************************/
    int numCompatible = 0;
                
    std::vector<PointAndValue<unsigned int> >::const_iterator firstEndpointIter;
    std::vector<PointAndValue<unsigned int> >::const_iterator secondEndpointIter;
                
    for (firstEndpointIter = firstEndpoint.myTree->getMyData()->begin();
         firstEndpointIter != firstEndpoint.myTree->getMyData()->end();
         firstEndpointIter++) {
                    
        for (secondEndpointIter = secondEndpoint.myTree->getMyData()->begin();
             secondEndpointIter != secondEndpoint.myTree->getMyData()->end();
             secondEndpointIter++) {

            // jmyers - this should be identical to some of the code in 
            // buildTracksAddToResults
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
    workItemFile myItem;
    myItem.endpoints.clear();
                
    //first item is first endpoint
    unsigned int treeId = firstEndpoint.myTime.getImageId();
    unsigned int nodeId = firstEndpoint.myTree->getId();
    myItem.endpoints.push_back(std::make_pair(treeId, nodeId));

    //second item is second endpoint
    treeId = secondEndpoint.myTime.getImageId();
    nodeId = secondEndpoint.myTree->getId();
    myItem.endpoints.push_back(std::make_pair(treeId, nodeId));

    //add all endpoints
    for(unsigned int m=0; m < newSupportNodes.size(); ++m){
        treeId = newSupportNodes.at(m).myTime.getImageId();
        nodeId = newSupportNodes.at(m).myTree->getId();
        myItem.endpoints.push_back(std::make_pair(treeId, nodeId));
    }

    //create workItemFile struct to go on list
    myItem.numNodes = 2 + newSupportNodes.size();
    myItem.numCompatibleEndpoints = numCompatible;
    myItem.timeUnits = myItem.numNodes * myItem.numCompatibleEndpoints;
    totalTimeUnits += myItem.timeUnits;

    assignment.at(globalNextWorker).push_back(myItem);

    ++globalNextWorker;
    if(globalNextWorker >= (numProcs-1)){
        globalNextWorker = 0;
    }

    ++totalWorkItems;

    //distribute workload periodically
    if(totalWorkItems >= MAX_WORK_ITEMS){

        distributeCurrentWorkload(assignment);/*, allDetections, allTracklets, searchConfig, 
                                                trackletTimeToTreeMap,
                                                totalTimeUnits);*/

        //update and reset assignment-specific variables
        totalWorkItems = 0;
        totalTimeUnits = 0;

        //clear our our work item cache and start afresh
        for(unsigned int i=0; i < assignment.size(); ++i){
            assignment.at(i).clear();
        }
        //assignment.clear(); //MATT NEW ADD 12/16/10

        //tell the workers that this round of work item distribution is complete and
        //collect any tracks that were generated
        killSentinel ta;
        ta.numProcs = numProcs;
        ta.sentinel = -2;
        killProcs(ta, allDetections);

        //only do MAX_GENERATIONS generations for measuring purposes
        std::cerr << "Generation " << numGenerations << std::endl;
        if(numGenerations >= MAX_GENERATIONS){
            return;
        }
        ++numGenerations;
    }

}



/***************************************************************************
 ***************************************************************************
 **   END DISTRIBUTED LINKTRACKLETS FUNCTION DEFINITIONS
 ***************************************************************************
 ***************************************************************************/
















/* **********************************************************************
* LINKTRACKLETS: The actual algorithm implementation *
* all this code by jmyers and should be identical to that 
* in trunk linkTracklets.cc
* **********************************************************************/


/*
 * update position and velocity given acceleration over time.
 */
void modifyWithAcceleration(double &position, double &velocity, 
                            double acceleration, double time)
{
    // use good ol' displacement = vt + .5a(t^2) 
    double newPosition = position + velocity*time 
        + .5*acceleration*(time*time);
    double newVelocity = velocity + acceleration*time;
    position = newPosition;
    velocity = newVelocity;
}




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










/***************************************************************************
 *                                                                         *
 *   LINKTRACKLETS FUNCTIONS WITH DISTRIBUTION LOGIC                       *
 *                                                                         *
 *   These segments a mixture of code by jmyers and mgcleveland.           *
 *   They should be functionally equivalent to the non-distributed         *
 *   versions but they will certainly not be identical.                    *
 *                                                                         *
 ***************************************************************************/




/*
 * this is called when all endpoint nodes (i.e. model nodes) and
 * support nodes are leaves.  model nodes and support nodes are
 * expected to be mutually compatible.
 */
/* the original was a void but it appears matt is returning some value
 * required for distribution logic. */
unsigned int buildTracksAddToResults(
    const std::vector<MopsDetection> &allDetections,
    const std::vector<Tracklet> &allTracklets,
    linkTrackletsConfig searchConfig,
    TreeNodeAndTime &firstEndpoint,
    TreeNodeAndTime &secondEndpoint,
    const std::vector<treeIdNodeIdPair> &supportNodes,
    TrackSet & results,

    // the following are arguments needed for distribution work.
    int &cacheSize, cachedTreeNode *nodeCache,
    std::map<ImageTime, TrackletTree > &trackletTimeToTreeMap,
    const std::vector<std::vector<treeIdNodeIdPair> > &finalEndpointOrder,
    unsigned long long &pageFaults, unsigned int &myNumCompatible,
    unsigned int setNum)
{

    /* matt cleveland distribution setup */
    //load support nodes from cache
    std::vector<std::vector<unsigned int> > supportNodesFromCache;
    unsigned int currPf = 0;
    unsigned int count =0;

    unsigned int nodeNum = 3; //start as first node offset

    std::vector<treeIdNodeIdPair>::const_iterator supportNodeIter;

    for (supportNodeIter = supportNodes.begin(); supportNodeIter != supportNodes.end();
         supportNodeIter++) {
        std::vector<unsigned int> dummy;
        supportNodesFromCache.push_back(dummy);
        treeIdNodeIdPair tini = *supportNodeIter;
        int i = loadNodeFromCache(tini.treeId, tini.nodeId, cacheSize, 
                                  nodeCache, trackletTimeToTreeMap,
                                  finalEndpointOrder, pageFaults, currPf,
                                  setNum, nodeNum);
        const std::vector<PointAndValue <unsigned int> > * curSupportNodeData;
        if(!nodeCache[i].node.myTree->isLeaf()){
            throw LSST_EXCEPT(BadParameterException,
			std::string(__FUNCTION__) + 
                              std::string(": received non-leaf node as support node."));
        }
        
        std::vector<PointAndValue <unsigned int> >::const_iterator supportPointIter;
        curSupportNodeData = nodeCache[i].node.myTree->getMyData(); 
        for (supportPointIter  = curSupportNodeData->begin(); 
             supportPointIter != curSupportNodeData->end();
             supportPointIter++) {
            //candidateTrackletIDs.push_back(supportPointIter->getValue());
            supportNodesFromCache.at(count).push_back(supportPointIter->getValue());
        }
        ++count;
    }
    
    /* begin standard jmyers logic like in trunk */
    std::vector<PointAndValue<uint> >::const_iterator firstEndpointIter;
    std::vector<PointAndValue<uint> >::const_iterator secondEndpointIter;
    std::vector<PointAndValue <uint> >::const_iterator supportPointIter;

    for (firstEndpointIter = firstEndpoint.myTree->getMyData()->begin();
         firstEndpointIter != firstEndpoint.myTree->getMyData()->end();
         firstEndpointIter++) {

        for (secondEndpointIter = secondEndpoint.myTree->getMyData()->begin();
             secondEndpointIter != secondEndpoint.myTree->getMyData()->end();
             secondEndpointIter++) {

            /* 
             * figure out the rough quadratic track fitting the two endpoints.
             * if error is too large, quit. Otherwise, choose support points
             * from the support nodes, using best-fit first, and ignoring those
             * too far off the line.  If we get enough points, return a track.
             *
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
	      // vars for treIdNodeId ==> TrackletTreeNode * conversion
	      // mgc 3/26/11
	      lsst::mops::TrackletTree     * trackletTree;
	      lsst::mops::TrackletTreeNode * treeNode;
	      // end mgc

                std::vector<uint> candidateTrackletIds;
                // put all support tracklet Ids in curSupportNodeData,
                // then call addDetectionsCloseToPredictedPositions
                for (supportNodeIter = supportNodes.begin(); 
                     supportNodeIter != supportNodes.end();
                     supportNodeIter++) {
                    const std::vector<PointAndValue <uint> > * 
                        curSupportNodeData;
		    
		    // mgc 3/26/11 -- get the tracklettree/node for this treeIdNodeId
		    trackletTree = findTreeById(supportNodeIter->treeId, trackletTimeToTreeMap);
		    
		    if (trackletTree != NULL) {
		      treeNode  = getNodeByIDAndTime(trackletTree, supportNodeIter->nodeId);

		      if (treeNode != NULL) {
		      
			if (!treeNode->isLeaf()) { // end mgc
			  throw LSST_EXCEPT(BadParameterException,
					    std::string(__FUNCTION__) + 
					    std::string(
							": received non-leaf node as support node."));
			}
			//curSupportNodeData = supportNodeIter->myTree->getMyData(); 
			curSupportNodeData = treeNode->getMyData(); //mgc
			for (supportPointIter  = curSupportNodeData->begin(); 
			     supportPointIter != curSupportNodeData->end();
			     supportPointIter++) {
			  
			  candidateTrackletIds.push_back(
							 supportPointIter->getValue());
			}
		      }
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

            }
        }    
    }

  
    /* this return added by mgcleveland (why?) */
    return currPf;
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
                      linkTrackletsConfig searchConfig,
                      TreeNodeAndTime &firstEndpoint,
                      TreeNodeAndTime &secondEndpoint,
                      std::vector<TreeNodeAndTime> &supportNodes,
                      double accMinRa, double accMaxRa, 
                      double accMinDec, double accMaxDec,
                      int iterationsTillSplit,
                      /* these items added by matt cleveland to allow distribution */
                      int numProcs, 
                      std::map<ImageTime, TrackletTree > &trackletTimeToTreeMap,
                      std::vector<std::vector<workItemFile> > &assignment)
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
                
                /**************************************************************
                 **************************************************************
                 **  DISTRIBUTED VERSION DOES NOT CALL buildTracksAddTResults
                 **  DIRECTLY, IT CREATES A SET OF WORK ITEMS AND THEN 
                 **  DISTRIBUTES THEM AFTER ANALYSIS
                 **************************************************************
                 **************************************************************/

                addWorkItem(allDetections, allTracklets, searchConfig,
                            firstEndpoint, secondEndpoint, newSupportNodes,
                            numProcs, trackletTimeToTreeMap, assignment);
                
            }
            else  {

                // jmyers: all this remaining code should be identical
                // to the sequential version except for the three
                // added arguments to each doLinkingRecurse call.

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
                                         iterationsTillSplit,
                                         numProcs,
                                         trackletTimeToTreeMap,
                                         assignment);
                        //std::cout << "Returned from recursion on
                        //left child of first endpoint.\n";
                    }
                    
                    if (firstEndpoint.myTree->hasRightChild())
                    {
                        TreeNodeAndTime newTAT(
                            firstEndpoint.myTree->getRightChild(), 
                            firstEndpoint.myTime);
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
                                         iterationsTillSplit,
                                         numProcs,
                                         trackletTimeToTreeMap,
                                         assignment);
					 
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
                                         iterationsTillSplit,
                                         numProcs,
                                         trackletTimeToTreeMap,
                                         assignment);
                        //std::cout << "Returned from recursion on
                        //left child of second endpoint.\n";
                    }
                    
                    if (secondEndpoint.myTree->hasRightChild())
                    {
                        TreeNodeAndTime newTAT(
                            secondEndpoint.myTree->getRightChild(), 
                            secondEndpoint.myTime);
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
                                         iterationsTillSplit,
                                         numProcs,
                                         trackletTimeToTreeMap,
                                         assignment);
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
               linkTrackletsConfig searchConfig,
               std::map<ImageTime, TrackletTree > &trackletTimeToTreeMap,
               // the following arguments added by matt for distribution work
	       int numProcs,
	       std::vector<std::vector<workItemFile> > &assignment)
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
                                         ITERATIONS_PER_SPLIT,
                                         numProcs, 
                                         trackletTimeToTreeMap,
                                         assignment);

                        if (searchConfig.myVerbosity.printStatus) {
                            time_t rawtime;
                            
                            struct tm * timeinfo;
                            std::cout << "That iteration took " 
                                      << timeSince(iterationTime) << " seconds. "
                                      << std::endl;
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
                        const linkTrackletsConfig &searchConfig,
                        uint numProcs) {

  std::string resultsDir("results");
  //boost::filesystem::remove_all(resultsDir);
  // make the directory in which worker procesors write their result Tracks
  mkdir(resultsDir.c_str(), 0777);

    /* system setup for distributed linkTracklets */
    srand(time(NULL));
    //work item set assignment creation
    std::vector<std::vector<workItemFile> > assignment;
    std::vector<workItemFile> dummyVector;
    for(uint i=1; i < numProcs; ++i){
        assignment.push_back(dummyVector);
    }
    /* end distributed linkTracklets setup */
    
    //TrackSet * toRet;
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
        std::cout << "Making detections contiguous.\n";
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

    doLinking(allDetections, 
              queryTracklets, 
              searchConfig, 
              trackletTimeToTreeMap, 
              numProcs,
              assignment);
    if (searchConfig.myVerbosity.printStatus) {
        std::cout << "Finished linking.\n";
    }    
  



    /* distributed linkTracklets cleanup */
    //distribute any remaining work items
    if(totalWorkItems > 0){
        std::cerr << timestamp() << "Distributing remaining work items." << std::endl;
        distributeCurrentWorkload(assignment);
    }
    else{
        std::cerr << timestamp() << 
            "There are no more work items for this input set." << std::endl;
    }
    
    std::cerr << timestamp() << "Killing final procs" << std::endl;
    killSentinel ta;
    ta.numProcs = numProcs;
    ta.sentinel = -1;
    killProcs(ta, allDetections);
    std::cerr << timestamp() << "Final procs killed" << std::endl;
    
    //wait for all threads to clear
    unsigned int SLEEP_SECS = 10;
    while(!stopAnnealing){
        std::cerr << "Dont stop annealing, hold on for " << 
            SLEEP_SECS << " seconds." << std::endl;
        sleep(SLEEP_SECS);
    }
    std::cerr << timestamp() << "Master returning with toRet size " << toRet->size() << std::endl;
    /* end cleanup stuff for distributed linkTracklets */

    return toRet;
}

}} /* close namespaces */


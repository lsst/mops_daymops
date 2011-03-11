// -*- LSST-C++ -*-
/* jonathan myers */

#ifndef LSST_LINKTRACKLETS_H_
#define LSST_LINKTRACKLETS_H_

#include <vector>
#include <time.h>
#include <map>

//#include "../KDTree.h"
//#include "
#include "../../include/lsst/mops/MopsDetection.h"
#include "../../include/lsst/mops/Tracklet.h"
#include "../../include/lsst/mops/KDTree.h"
#include "../../include/lsst/mops/KDTreeNode.h"
//#include "../Tracklet.h"
#include "../../include/lsst/mops/Track.h"
#include "../../include/lsst/mops/TrackSet.h"


/*
 * linkTrackletsConfig is a simple class which holds the many configurable
 * parameters for running linkTracklets.  The default values should be 
 * fine for LSST.
 */


class linkTrackletsConfig {
public:

    linkTrackletsConfig() 
        {   maxRAAccel = .02; 
            maxDecAccel = .02; 
            detectionLocationErrorThresh = .002; 
            minEndpointTimeSeparation = 2; 
            minSupportToEndpointTimeSeparation = .5;
            minSupportTracklets = 1;
            quadraticFitErrorThresh = 0.;
            minDetectionsPerTrack = 6;

            velocityErrorThresh = .00025;//.5;  per Jon's email 5/20/10
	    //.5 turned out to miss some potential tracks;  11/23/10 - setting back to .00025 to test out performance
        }

    /* acceleration terms are in degrees/(day^2)
       
       objects which accelerate faster than these maximum acceleration params
       are not reported (or considered)
     */
    double maxRAAccel;
    double maxDecAccel;


    /* detection error thresh is the upper bound on observational error for
       Detections; this has repercussions for what Detections are included in a
       returned track as well as the behavior of the searching itself (we will
       prune off regions of space based on their location, so error threshold
       will factor in here.)
     */
    double detectionLocationErrorThresh;


    /* quadratic fit error thresh: 

       this is used to determine how much error we will accept from the
       quadratic fitting process.  If a detection is within
       detectionLocationErrorThresh + quadraticFitErrorThresh of the predicted
       location, it will be considered for addition to a track
     */
    double quadraticFitErrorThresh;


    /*
      minSupportToEndpointTimeSeparation: the minimum time between endpoints and
      support points.  normally, orbit fits need several days of observations in
      order to be useful, so we will ignore support points which happen right
      after the actual detection. This also speeds up the search a little.
    */
    double minSupportToEndpointTimeSeparation;
    /*
      min time between the first and last point of a tracklet in order to be
      considered for a track
     */
    double minEndpointTimeSeparation;
    /*
      minimum number of *support* tracklets required for a track.  note that the
      track will have two "endpoint" tracklets PLUS this many support tracklets.
     */
    double minSupportTracklets;
    
    /*
      the minimum number of detections per track
     */
    unsigned int minDetectionsPerTrack;

    /*
      velocity error thresh: this is used when calculating whether two regions
      could contain the same object.  this is the maximum believable discrepancy
      between *actual* velocity and best-fit velocity, measured in deg/day.
     */
    double velocityErrorThresh;

};





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
  TreeNodeAndTime(lsst::mops::KDTreeNode<unsigned int> * tree, ImageTime i) {
    myTree = tree;
    myTime = i;
  }
  lsst::mops::KDTreeNode <unsigned int> * myTree;
  ImageTime myTime;
};


/***************************************************************************
 * This is the data that is stored in the node cache.  A node cache entry
 * contains the unique identifiers for the TreeNodeAndTime object (treeId,
 * nodeId) and the TreeNodeAndTime object itself.
 ***************************************************************************/
typedef struct cn{
  unsigned int treeId, nodeId;
  TreeNodeAndTime node;
  bool inUse;
  time_t timeStamp;
} cachedNode;

typedef struct ta{
  int numProcs;
  int sentinel;
}threadArgs;

/**************************************************
 * This struct is used for keeping track of the
 * work item and meta data about them, which
 * will be used when the work is distributed.
 **************************************************/
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


/* queryTracklets are non-const because we set their velocityRA and velocityDec fields. 
   otherwise queryTracklets will not be changed. */
void linkTracklets(const std::vector<lsst::mops::MopsDetection> &allDetections,
		   std::vector<lsst::mops::Tracklet> &queryTracklets,
		   linkTrackletsConfig searchConfig, int numProcs,
		   lsst::mops::TrackSet &tracks);

void waitForTask(int rank,
		 const std::vector<lsst::mops::MopsDetection> &allDetections, //from MAIN
		 std::vector<lsst::mops::Tracklet> &allTracklets, //from MAIN
		 linkTrackletsConfig searchConfig  /*from MAIN*/);
	

//only declared for unit testing
void getBestFitVelocityAndAcceleration(std::vector<double> positions, const std::vector<double>&times,
                                       double & velocity, double &acceleration, double &position0);

// note that position and velocity will be MODIFIED. 
void modifyWithAcceleration(double &position, double &velocity, 
                            double acceleration, double time);


/***************************************************************************
 * DLT function prototypes
 ***************************************************************************/
void 
distributeCurrentWorkload(std::vector<std::vector<workItemFile> > &assignment);

void 
*killProcs(threadArgs);

int
leastRecentlyUsed(cachedNode *nodeCache, int cacheSize);

int 
findFarthestNodeIndex(cachedNode *nodeCache, int cacheSize,
		      const std::vector<std::vector<treeIdNodeIdPair> > &finalEndpointOrder,
		      unsigned int setStart, unsigned int workItemStart);

lsst::mops::KDTree<unsigned int> *
findTreeById(const unsigned int id, 
	     std::map<ImageTime, lsst::mops::KDTree<unsigned int> > &treeMap);

ImageTime 
findImageTimeForTree(const unsigned int id, 
			       const std::map<ImageTime, lsst::mops::KDTree<unsigned int> > &treeMap);

lsst::mops::KDTreeNode<unsigned int> *
getNodeByIDAndTime(lsst::mops::KDTree<unsigned int> *myTree, unsigned int nodeId);

int 
loadNodeFromCache(unsigned int treeId, unsigned int nodeId, int &cacheSize, cachedNode *nodeCache,
		  std::map<ImageTime, lsst::mops::KDTree <unsigned int> > &trackletTimeToTreeMap,
		  const std::vector<std::vector<treeIdNodeIdPair> > &finalEndpointOrder,
		  unsigned long long &pageFaults);

#endif 

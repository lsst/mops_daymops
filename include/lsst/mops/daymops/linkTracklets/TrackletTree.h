// -*- LSST-C++ -*-



/* jmyers 2/10/11
 *
 * Much like KDTree, TrackletTree is the "head" only, and TrackletTreeNode does
 * the real work. See TrackletTreeNode.h for more comments.
 *
 * Tracklet trees are for holding tracklet *rooted in a single image* and are
 * partitioned by RA, Dec, RAv, Decv.  The nodes have potentially *overlapping*
 * bounds, because the tracklets have error bounds on them.
 * 
 * This is very similar to KDTree, so check KDTree.h for comments. Unlike
 * KDTree, TrackletTrees are for use with linkTracklets - they have potentially
 * overlapping nodes, which are extended by error bounds in the tracklets.
 */


#ifndef TRACKLET_TREE_H
#define TRACKLET_TREE_H

#include <vector>
#include <cmath>
#include <algorithm>

#include "lsst/mops/common.h"  
#include "lsst/mops/Exceptions.h"
#include "lsst/mops/PointAndValue.h"
#include "lsst/mops/KDTree.h"
#include "lsst/mops/daymops/linkTracklets/TrackletTreeNode.h"

#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/TrackletVector.h"



namespace lsst {
namespace mops {

    
    class TrackletTree: public KDTree<unsigned int> {
    public:

        /* default constructor - don't instantiate with any data.         
         */
        TrackletTree();

        /*
         * Give the MopsDetections and Tracklets rooted in a given image.
         * Builds a TrackletTree on (RA, Dec, RAv, Decv). 
         *
         * The bounds of the tree nodes are extended by the positional error
         * bars given, and the max velocities are extended as well using
         * delta_time and positional error.
         *
         * maxLeafSize should be a positive integer.
         */
        TrackletTree(const std::vector<MopsDetection> &allDetections,
                     const TrackletVector &thisTreeTracklets,
                     double positionalErrorRa, double positionalErrorDec,
                     unsigned int maxLeafSize);

        /* copy constructor */
        TrackletTree(const TrackletTree &source);

        /* 
         * populates the tree with given data.  See second constructor.
         */
        void buildFromData(const std::vector<MopsDetection> &allDetections,
                           const TrackletVector &thisTreeTracklets,
                           double positionalErrorRa, double positionalErrorDec,
                           unsigned int maxLeafSize);


        // TBD: so... since KDTree has a KDTreeNode<T> myRoot, what happens
        // here? Do normal KDTree operations which use root simply fail
        // spectacularly?  probably, yes.

        TrackletTreeNode *myRoot;

    };






TrackletTree::TrackletTree(const std::vector<MopsDetection> &allDetections,
                           const TrackletVector &thisTreeTracklets,
                           double positionalErrorRa, double positionalErrorDec,
                           unsigned int maxLeafSize)
{
    setUpEmptyTree();
    buildFromData(allDetections, thisTreeTracklets, positionalErrorRa,
                  positionalErrorDec, maxLeafSize);

}




void TrackletTree::buildFromData(const std::vector<MopsDetection> &allDetections,
                                 const TrackletVector &thisTreeTracklets,
                                 double positionalErrorRa, double positionalErrorDec,
                                 unsigned int maxLeafSize)
{
    // need to set up fields used by KDTree just the way KDTree would; then 
    // create our set of child TrackletTreesNodes
    myK = 4;

    if (thisTreeTracklets.size() > 0) 
    {

        std::vector<PointAndValue <unsigned int> > parameterizedTracklets;
        std::vector<double> pointsUBounds, pointsLBounds;

        if (maxLeafSize < 1) {
            throw LSST_EXCEPT(BadParameterException, 
                              "EE: KDTree: max leaf size must be strictly positive!\n");
        }


        // need to convert tracklets to a parameterized format (RA_0, Dec_0, RAv, Dec_v)
        // calculate UBounds, LBounds while we're at it

        // TBD:

        
        // build root TrackletTreeNode

        /* create the root of the tree (and the rest of the tree
         * recursively), save it to private var. */
        unsigned int idCounter = 0;
        myRoot = new TrackletTreeNode(parameterizedTracklets, 
                                      positionalErrorRa, positionalErrorDec, maxLeafSize, 0,
                                      pointsUBounds, pointsLBounds, 
                                      idCounter);
        // don't set hasData until now, when the tree is actually built.
        mySize = idCounter;

        // *read* UBounds, UBounds from TrackletTreeNode (They have been
        // *updated) and set our own

  
        hasData = true;
    }
}













}} // close namespace lsst::mops

#endif

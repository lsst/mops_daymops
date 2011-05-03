// -*- LSST-C++ -*-


#ifndef TRACKLET_TREE_NODE_H
#define TRACKLET_TREE_NODE_H

#include <iostream>

#include "lsst/mops/BaseKDTreeNode.h"
#include "lsst/mops/Exceptions.h"
#include "lsst/mops/LinkageVector.h"


/***************************************************************************
 *
 * trackletTrees are 4D KDTrees which are for holding tracklets, and not
 * anything else.  Unlike normal KDTreeNodes, they are to be used by
 * linkTracklets/ the variable trees algorithm (Kubica), and they are NOT to be
 * used for range searching.  Each node is extended by the max positional /
 * velocity error of its children, so trackletTree nodes can overlap!
 * 
 ****************************************************************************/

namespace lsst {
namespace mops {


    class TrackletTreeNode: 
        public BaseKDTreeNode<unsigned int, TrackletTreeNode> {
    public: 
        friend class BaseKDTreeNode<unsigned int, TrackletTreeNode>;

        /* linkages should be a series of pointAndValues and should
         * hold a bunch of elements like:
         * 
         * (RA, Dec, RAv, DecV, deltaTime)
         * 
         * or
         *
         * (RA, Dec, RAv, DecV, RAAcc, DecAcc, deltaTime)
         *
         * of the tracklet, mapped to the ID of the tracklet.
         * 
         * Max allowable positional error is specified by the user;
         * this is used along with tracklet delta-time to calculate
         * the velocity error.
         *
         * if useMedian = false, splitWidest=true and widths is nonempty,
         * we split up values like C linkTracklets does. 
         */


        TrackletTreeNode(
            const std::vector<MopsDetection> allDetections,
            const std::vector<PointAndValue <unsigned int> > &linkages, 
            LinkageVector* allLinkages,
            double positionalErrorRa, 
            double positionalErrorDec,
            unsigned int maxLeafSize, 
            unsigned int myAxisToSplit, 
            const std::vector<double> &widths,
            unsigned int &lastId,
            bool useMedian=false,
            bool splitWidest=true);
        
    
        const unsigned int getNumVisits() const;
        void addVisit();


        // these are to be used by linkTracklets.
        bool hasLeftChild() const;
        bool hasRightChild() const;
        TrackletTreeNode * getLeftChild();
        TrackletTreeNode * getRightChild();        
        // return IDs of objects with true tracks at or below this node
        const std::set<int> * getMyObjects();

        /* getMyData returns PointsAndValues; the Values are indexes
         * into a vector. getDataParentVec is a pointer to that
         * vector.  */
        LinkageVector* getDataParentVec() const;

        const std::vector<PointAndValue <unsigned int> > * getMyData() const;
        bool isLeaf() const;




        // return true iff this node OR ITS CHILDREN holds the
        // tracklet w. ID t. This is crazy slow on large trees (due to
        // the fact that it does an undirected depth-first
        // traversal) but helpful for debugging.
        bool hasTracklet(unsigned int t);

    protected:

        LinkageVector* dataParentVec;

        /* after the "real" constructor is called at the root, and the
         * BaseKDTree constructor sets up the children, this function
         * is used to do a post-traversal of the children and update
         * their error bounds.
         */
        void recalculateBoundsWithError(double positionalErrorRa,
                                        double positionalErrorDec);

    private:
        unsigned int numVisits;
        std::set<int> myObjects;
    };



}} // close namespace lsst::mops

#endif

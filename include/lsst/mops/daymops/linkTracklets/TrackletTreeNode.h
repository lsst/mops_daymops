// -*- LSST-C++ -*-


#ifndef TRACKLET_TREE_NODE_H
#define TRACKLET_TREE_NODE_H

#include <iostream>

#include "lsst/mops/BaseKDTreeNode.h"
#include "lsst/mops/Exceptions.h"


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


    class TrackletTreeNode: public BaseKDTreeNode<unsigned int, TrackletTreeNode> {
    public: 
        friend class BaseKDTreeNode<unsigned int, TrackletTreeNode>;

        /* tracklets is a series of pointAndValues and should hold a
         * bunch of elements like:
         * 
         * (RA, Dec, RAv, DecV, deltaTime)
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
            const std::vector<PointAndValue <unsigned int> > &tracklets, 
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

        // return true iff this node OR ITS CHILDREN holds the tracklet t
        bool hasTracklet(unsigned int t);

        // these are to be used by linkTracklets.
        bool hasLeftChild() const;
        bool hasRightChild() const;
        TrackletTreeNode * getLeftChild();
        TrackletTreeNode * getRightChild();

        const std::vector<PointAndValue <unsigned int> > * getMyData() const;
        bool isLeaf() const;


    protected:

        /* after the "real" constructor is called at the root, and the
         * BaseKDTree constructor sets up the children, this function
         * is used to do a post-traversal of the children and update
         * their error bounds.
         */
        void recalculateBoundsWithError(double positionalErrorRa,
                                        double positionalErrorDec);

    private:
        unsigned int numVisits;
    };



}} // close namespace lsst::mops

#endif

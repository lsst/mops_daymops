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
         */


        TrackletTreeNode(
            const std::vector<PointAndValue <unsigned int> > &tracklets, 
            double positionalErrorRa, 
            double positionalErrorDec,
            unsigned int maxLeafSize, 
            unsigned int myAxisToSplit, 
            const std::vector<double> &UBounds,
            const std::vector<double> &LBounds, unsigned int &lastId);
        
    
        const unsigned int getNumVisits() const;
        void addVisit();
        

        // these are to be used by linkTracklets.
        bool hasLeftChild() const;
        bool hasRightChild() const;
        TrackletTreeNode * getLeftChild();
        TrackletTreeNode * getRightChild();

        const std::vector<PointAndValue <unsigned int> > * getMyData() const;
        bool isLeaf() const;


    protected:
        /* Our "real" constructor calls the BaseKDTreeNode constructor
         * then does some followup.  The BaseKDTreeNode constructor
         * needs to call *our* constructor with this form (which in
         * turn just calls the BaseKDTree constructor of the same
         * form)!
         */
        TrackletTreeNode(
            const std::vector<PointAndValue <unsigned int> > &pointsAndValue,
            unsigned int k, 
            unsigned int maxLeafSize, 
            unsigned int myAxisToSplit, 
            const std::vector<double> & UBounds,
            const std::vector<double> & LBounds, 
            unsigned int &lastId);


    private:

        unsigned int numVisits;
    };



}} // close namespace lsst::mops

#endif

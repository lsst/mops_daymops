// -*- LSST-C++ -*-


#ifndef TRACKLET_TREE_NODE_H
#define TRACKLET_TREE_NODE_H

#include <iostream>

#include "lsst/mops/KDTree.h"
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


    class TrackletTreeNode: public KDTreeNode<unsigned int> {
    public: 

        /* tracklets is a series of pointAndValues and should hold a bunch of
         * elements like:
         * 
         * (RA, Dec, RAv, DecV, deltaTime)
         *
         * of the tracklet, mapped to the ID of the tracklet.  
         * 
         * Max allowable positional error is specified by the user; this is used
         * along with tracklet delta-time to calculate the velocity error.
         */


        TrackletTreeNode(std::vector<PointAndValue <unsigned int> > tracklets, 
                         double positionalErrorRa, double positionalErrorDec,
                         unsigned int maxLeafSize, 
                         unsigned int myAxisToSplit, std::vector<double> UBounds,
                         std::vector<double>LBounds, unsigned int &lastId);
        
        /*
         * this addReference function is to be used by KDTree *ONLY*.
         * This way multiple KDTrees can share the same nodes without
         * expensive copying (and since trees are static, this is fine.)
         *
         * it is up to KDTree to handle the creation and deletion of these nodes
         * based on refcount.
         */
        // TBD: make KDTreeNode a friend class and make these protected
        void addReference();

        void removeReference();

        unsigned int getRefCount();

        // WARNING: These search functions may not work! I'm including them
        // because C++ requires we have all the methods from the parent class.
        std::vector<PointAndValue <unsigned int> > 
        rangeSearch(std::vector<double> queryPt, 
                    double queryRange) const; 
    
        std::vector<PointAndValue <unsigned int> > 
        hyperRectangleSearch(const std::vector<double> &queryPt, 
                             const std::vector<double> &tolerances, 
                             const std::vector<GeometryType> &spaceTypesByDimension) const;


        // these are to be used by linkTracklets.
        bool hasLeftChild() const;
        bool hasRightChild() const;
        TrackletTreeNode * getLeftChild();
        TrackletTreeNode * getRightChild();
        /*return a pointer to a const vector of the per-axis upper bounds of
          this tree node.  The format is identical to the double vector in the
          pointsAndValues used to create this node.
         */

        // note that in tracklet trees, node bounds MAY OVERLAP
        const std::vector<double> *getUBounds() const;
        // same, but with the lower bounds rather than upper bounds.
        const std::vector<double> *getLBounds() const;
        const std::vector<PointAndValue <unsigned int> > * getMyData() const;
        bool isLeaf() const;
        
        void debugPrint(int depth) const;
    
        std::vector <TrackletTreeNode> myChildren;
        
        const unsigned int getNumVisits() const;
        const unsigned int getId() const;
        void addVisit();

        /* for a collection of points and values, find the median value
           along the given axis and return it

           ASSUMES all points have size > axis 
        */
        double getMedianByAxis(std::vector<PointAndValue<unsigned int> > pointsAndValues,
                               unsigned int axis);

        double maxByAxis(std::vector<PointAndValue<unsigned int> > pointsAndValues,
                         unsigned int axis);
        

    private:
        // helper function for extending UBounds/LBounds 
        void extendBounds(std::vector<double> &myBounds, const std::vector<double> &childBounds,
                          bool areUBounds);

    };




    TrackletTreeNode::TrackletTreeNode(std::vector<PointAndValue <unsigned int> > tracklets, 
                                       double positionalErrorRa, double positionalErrorDec,
                                       unsigned int maxLeafSize, 
                                       unsigned int myAxisToSplit, std::vector<double> UBounds,
                                       std::vector<double>LBounds, unsigned int &lastId)
        : KDTreeNode<unsigned int>(tracklets, 4, maxLeafSize, myAxisToSplit, UBounds, LBounds, lastId)
    {

        // KDTreeNode has already set up our refcounts, built our children or
        // made us a leaf, etc.
        if (KDTreeNode<unsigned int>::isLeaf()) {
            // extend UBounds, LBounds based on shortest delta_time in data.
            if (myData.size() <= 0) {
                LSST_EXCEPT(ProgrammerErrorException, 
                            "Unexpected condition: tree node is leaf, but has no data.");
            }
            // find shortest delta_timme in tracklets we own
            double shortestChildDt = tracklets.at(0).getPoint().at(4);
            for (unsigned int i = 0; i < tracklets.size(); i++) {
                double thisDt = tracklets.at(i).getPoint().at(4);
                if (thisDt < shortestChildDt) {
                    shortestChildDt = thisDt;
                }
            }
            double velocityErrorRa  = 2.0 * positionalErrorRa  / shortestChildDt;
            double velocityErrorDec = 2.0 * positionalErrorDec / shortestChildDt;
            UBounds.at(0) += positionalErrorRa;
            LBounds.at(0) -= positionalErrorRa;
            UBounds.at(1) += positionalErrorDec;
            LBounds.at(1) -= positionalErrorDec;
            UBounds.at(2) += velocityErrorRa;
            LBounds.at(2) -= velocityErrorRa;
            UBounds.at(3) += velocityErrorDec;
            LBounds.at(3) -= velocityErrorDec;
        }
        else {
            // extend our UBounds, LBounds by child max UBounds, LBounds.
            if (KDTreeNode<unsigned int>::hasLeftChild()) {
                extendBounds(UBounds, *(KDTreeNode<unsigned int>::getLeftChild()->getUBounds()),
                             true);
                extendBounds(LBounds, *(KDTreeNode<unsigned int>::getLeftChild()->getLBounds()),
                             false);

            }
            if (KDTreeNode<unsigned int>::hasRightChild()) {
                extendBounds(UBounds, *(KDTreeNode<unsigned int>::getRightChild()->getUBounds()),
                             true);
                extendBounds(LBounds, *(KDTreeNode<unsigned int>::getRightChild()->getLBounds()),
                             false);

            }
        }
    }



void TrackletTreeNode::extendBounds(std::vector<double> &myBounds, const std::vector<double> &childBounds,
                  bool areUBounds) 
{
    if (myBounds.size() != childBounds.size()) {
        LSST_EXCEPT(BadParameterException, "Expect argument vectors of same length.");
    }
    for (unsigned int i = 0; i < myBounds.size(); i++) {
        if (areUBounds) {
            if (childBounds.at(i) > myBounds.at(i)) {
                myBounds.at(i) = childBounds.at(i);
            }
        }
        else {
            if (childBounds.at(i) < myBounds.at(i)) {
                myBounds.at(i) = childBounds.at(i);
            }

        }
    }
}


}} // close namespace lsst::mops

#endif

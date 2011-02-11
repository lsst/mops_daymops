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

        /* TBD: I haven't actually tested that these work, but a quick glance
         * over the code in KDTreeNode says they probably will.  No matter
         * either way; linkTracklets doesn't use them.
         */
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
        const std::vector<double> *getLBounds() const;

        const std::vector<PointAndValue <unsigned int> > * getMyData() const;
        bool isLeaf() const;
        
        void debugPrint(int depth) const;
    
        std::vector <TrackletTreeNode> myChildren;
        
        const unsigned int getNumVisits() const;
        const unsigned int getId() const;
        void addVisit();

        double getMedianByAxis(std::vector<PointAndValue<unsigned int> > pointsAndValues,
                               unsigned int axis);

        double maxByAxis(std::vector<PointAndValue<unsigned int> > pointsAndValues,
                         unsigned int axis);
        

    private:
        // helper function for extending UBounds/LBounds 
        void extendBounds(std::vector<double> &myBounds, 
                          const std::vector<double> &childBounds,
                          bool areUBounds);

    };




    TrackletTreeNode::TrackletTreeNode(std::vector<PointAndValue <unsigned int> > tracklets, 
                                       double positionalErrorRa, double positionalErrorDec,
                                       unsigned int maxLeafSize, 
                                       unsigned int myAxisToSplit, 
                                       std::vector<double> UBounds,
                                       std::vector<double>LBounds, unsigned int &lastId)

    /* build a 4-D tree using these parameterized tracklets, which are really
     * 5-dimensional. KDTreeNode will politely ignore the last item, which is
     * expected to be delta_time, and just partition on the 4 dimensions we care
     * about: RA, Dec, RAv, Decv.
     */

        : KDTreeNode<unsigned int>(tracklets, 4, maxLeafSize, 
                                   myAxisToSplit, UBounds, LBounds, lastId)
    {

        /* The call to the KDTreeNode constructor has already set up our
         * refcounts, built our children or made us a leaf, etc. Now just extend
         * UBounds, LBounds. */
        
        if (KDTreeNode<unsigned int>::isLeaf()) {

            // extend UBounds, LBounds with knowledge of velocity error.
            if (myData.size() <= 0) {
                LSST_EXCEPT(ProgrammerErrorException, 
                            "Unexpected condition: tree node is leaf, but has no data.");
            }

            // extend UBounds, LBounds in RA, Dec by position error.
            UBounds.at(0) += positionalErrorRa;
            LBounds.at(0) -= positionalErrorRa;
            UBounds.at(1) += positionalErrorDec;
            LBounds.at(1) -= positionalErrorDec;

            // find min/max RA, Dec velocities after accounting for error.

            for (unsigned int i = 0; i < myData.size(); i++) {
                std::vector<double> trackletPoint = myData.at(i).getPoint();
                double trackletRaV  = trackletPoint.at(2);
                double trackletDecV = trackletPoint.at(3);
                double thisDt = trackletPoint.at(4);

                double maxRaV = trackletRaV + 2.0 * positionalErrorRa / thisDt;
                double minRaV = trackletRaV - 2.0 * positionalErrorRa / thisDt;
                double maxDecV = trackletDecV + 2.0 * positionalErrorDec / thisDt;
                double minDecV = trackletDecV - 2.0 * positionalErrorDec / thisDt;
                
                if (UBounds.at(2) < maxRaV) {
                    UBounds.at(2) = maxRaV;
                }
                if (UBounds.at(3) < maxDecV) {
                    UBounds.at(3) = maxDecV;
                }
                if (LBounds.at(2) > minRaV) {
                    LBounds.at(2) = minRaV;
                }
                if (LBounds.at(3) < minDecV) {
                    LBounds.at(3) = minDecV;
                }
            }

        }
        else {
            // extend our UBounds, LBounds by child max UBounds, LBounds.
            if (KDTreeNode<unsigned int>::hasLeftChild()) {
                extendBounds(UBounds, 
                             *(KDTreeNode<unsigned int>::getLeftChild()->getUBounds()),
                             true);
                extendBounds(LBounds, 
                             *(KDTreeNode<unsigned int>::getLeftChild()->getLBounds()),
                             false);

            }
            if (KDTreeNode<unsigned int>::hasRightChild()) {
                extendBounds(UBounds, 
                             *(KDTreeNode<unsigned int>::getRightChild()->getUBounds()),
                             true);
                extendBounds(LBounds, 
                             *(KDTreeNode<unsigned int>::getRightChild()->getLBounds()),
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

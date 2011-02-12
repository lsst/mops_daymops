#include "lsst/mops/daymops/linkTracklets/TrackletTreeNode.h"

namespace lsst { namespace mops {

// these are to be used by linkTracklets.
const unsigned int TrackletTreeNode::getNumVisits() const
{
    return numVisits;
}


void TrackletTreeNode::addVisit() 
{
    numVisits++;
}




bool TrackletTreeNode::hasLeftChild() const
{
    if (myChildren.size() >= 1) {
        return true;
    }
    return false;
}


bool TrackletTreeNode::hasRightChild() const
{
    if (myChildren.size() >= 2) {
        return true;
    }
    return false;
}


TrackletTreeNode * TrackletTreeNode::getLeftChild()
{
    if (!hasLeftChild()) {
        return NULL;
    }
    return &(myChildren.at(0));
}


TrackletTreeNode * TrackletTreeNode::getRightChild()
{
    if (!hasRightChild()) {
        return NULL;
    }
    return &(myChildren.at(1));
}


const std::vector<double> *TrackletTreeNode::getUBounds() const
{
    return &myUBounds;
}


const std::vector<double> *TrackletTreeNode::getLBounds() const
{
    return &myLBounds;
}


const std::vector<PointAndValue <unsigned int> > * 
TrackletTreeNode::getMyData() const
{
    if (!isLeaf()) {
        LSST_EXCEPT(BadParameterException,
                    "KDTreeNode got request for data, but is not a leaf node.");
    }
    return &myData;
}


bool TrackletTreeNode::isLeaf() const
{
    if (myChildren.size() == 0) 
    {
        return true;
    }
    return false;
}
    




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

        : BaseKDTreeNode<unsigned int, TrackletTreeNode>(tracklets, 4, maxLeafSize, 
                                                         myAxisToSplit, UBounds, LBounds, 
                                                         lastId)
    {

        numVisits = 0;

        /* The call to the KDTreeNode constructor has already set up our
         * refcounts, built our children or made us a leaf, etc. Now just extend
         * UBounds, LBounds. */
        
        if (isLeaf()) {

            // extend UBounds, LBounds with knowledge of velocity error.
            if (myData.size() <= 0) {
                LSST_EXCEPT(ProgrammerErrorException, 
                            "Unexpected condition: tree node is leaf, but has no data.");
            }

            // extend UBounds, LBounds in RA, Dec by position error.
            myUBounds.at(0) += positionalErrorRa;
            myLBounds.at(0) -= positionalErrorRa;
            myUBounds.at(1) += positionalErrorDec;
            myLBounds.at(1) -= positionalErrorDec;

            // find min/max RA, Dec velocities after accounting for error.

            for (unsigned int i = 0; i < myData.size(); i++) {
                std::vector<double> trackletPoint = myData.at(i).getPoint();
                double trackletRaV  = trackletPoint.at(2);
                double trackletDecV = trackletPoint.at(3);
                double thisDt = trackletPoint.at(4);

		double maxVelocityErrRa =  2.0 * positionalErrorRa  / thisDt;
		double maxVelocityErrDec = 2.0 * positionalErrorDec / thisDt;

                double maxRaV = trackletRaV + maxVelocityErrRa;
                double minRaV = trackletRaV - maxVelocityErrRa;
                double maxDecV = trackletDecV + maxVelocityErrDec;
                double minDecV = trackletDecV - maxVelocityErrDec;
                
                if (myUBounds.at(2) < maxRaV) {
                    myUBounds.at(2) = maxRaV;
                }
                if (myUBounds.at(3) < maxDecV) {
                    myUBounds.at(3) = maxDecV;
                }
                if (myLBounds.at(2) > minRaV) {
                    myLBounds.at(2) = minRaV;
                }
                if (myLBounds.at(3) > minDecV) {
                    myLBounds.at(3) = minDecV;
                }
            }

        }
        else {
            // extend our UBounds, LBounds by child max UBounds, LBounds.
            if (hasLeftChild()) {

                extendBounds(myUBounds, 
                             *(getLeftChild()->getUBounds()),
                             true);
                extendBounds(myLBounds, 
                             *(getLeftChild()->getLBounds()),
                             false);

            }
            if (hasRightChild()) {

                extendBounds(myUBounds, 
                             *(getRightChild()->getUBounds()),
                             true);
                extendBounds(myLBounds, 
                             *(getRightChild()->getLBounds()),
                             false);

            }
        }
    }
     


TrackletTreeNode::TrackletTreeNode(
     std::vector<PointAndValue <unsigned int> > pointsAndValues, 
     unsigned int k, 
     unsigned int maxLeafSize, 
     unsigned int myAxisToSplit, 
     std::vector<double> UBounds,
     std::vector<double>LBounds, 
     unsigned int &lastId)
{
     BaseKDTreeNode<unsigned int, TrackletTreeNode>::
	  BaseKDTreeNode
	  ( pointsAndValues, k, maxLeafSize, myAxisToSplit,
	    UBounds, LBounds, lastId);
}
     


}}//close lsst::Mops

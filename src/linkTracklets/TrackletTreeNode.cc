// -*- LSST-C++ -*-
/* jmyers 2/10/11 */

#include "lsst/mops/daymops/linkTracklets/TrackletTreeNode.h"
#include "lsst/mops/Exceptions.h"
#include "lsst/mops/MopsDetection.h"

#define uint unsigned int

namespace lsst { namespace mops {





LinkageVector* TrackletTreeNode::getDataParentVec() const
{
    return dataParentVec;
}

        
TrackletTreeNode::TrackletTreeNode(
    const std::vector<MopsDetection> allDetections,
    const std::vector<PointAndValue <unsigned int> > &tracklets, 
    LinkageVector* allLinkages,
    double positionalErrorRa,
    double positionalErrorDec,
    unsigned int maxLeafSize, 
    unsigned int myAxisToSplit, 
    const std::vector<double> &widths,
    unsigned int &lastId,
    bool useMedian,
    bool splitWidest)
{
    myRefCount = 1;
    dataParentVec = allLinkages;

    /* we expect either ra, dec, raV, decV, dt or
       ra, dec, raV, decV, raAcc, decAcc, dt.
     */
    myK = tracklets.at(0).getPoint().size() - 1;
    if ((myK != 4) && (myK != 6)) {
        throw LSST_EXCEPT(BadParameterException,
    "Don't know what to do with data of this dimensionality!\n");
    }

    lastId++;
    id = lastId;

    std::vector<double> rightChildUBounds, rightChildLBounds,     
        leftChildUBounds, leftChildLBounds;
    std::vector<PointAndValue <uint> > leftPointsAndValues, rightPointsAndValues;

    
    myUBounds.resize(myK);
    myLBounds.resize(myK);

    // need to calculate initial UBounds, LBounds for our data.
    for (uint i = 0; i < tracklets.size(); i++) {
        if (tracklets.at(i).getPoint().size() != myK + 1) {
            LSST_EXCEPT(ProgrammerErrorException, 
       "expected all tracklet/track points to be 5d or 7d: Ra, Dec, RaV, DecV, dt or Ra, Dec, RaV, DecV, RaAcc, DecAcc, dt.\n");
        }
        for (uint axis = 0; axis < myK; axis++) {
            double val = tracklets.at(i).getPoint().at(axis);
            if ((i == 0) || (val > myUBounds[axis])) {
                myUBounds[axis] = val;
            }
            if ((i ==0) || (val < myLBounds[axis])) {
                myLBounds[axis] = val;
            }
        }
    }
    

    if (tracklets.size() <= maxLeafSize) {
       // leaf case is easy.
        myData = tracklets;
        for (unsigned int i = 0; i != myData.size(); i++) {
            int obj= allLinkages->at(
                myData.at(i).getValue())
                ->getObjectId(allDetections);
            if (obj != -1) {
                myObjects.insert(obj);
            }
        }
    }

    else {
        // non-leaf case; we have much to do.
        double pivot;
        unsigned int nextAxis;


        // if not splitWidest, our parent gave us which axis to split.
        if (splitWidest) {
            // choose an axis to split, overwrite myAxisToSplit
            double maxWidth = -1;
            uint widestAxis = -1; 
            for (uint i = 0; i < myK; i++) {
                double width = (myUBounds[i] - myLBounds[i]) / widths[i];
                if (width < 0) {
                    LSST_EXCEPT(ProgrammerErrorException,
             " Got impossible width when building trackletTree\n");
                }
                if (width > maxWidth) {
                    maxWidth = width;
                    widestAxis = i;
                }
            }
            myAxisToSplit = widestAxis;
        }
        
        // split up data in our axis.
        if (useMedian) {
            // use the median.
            pivot = getMedianByAxis(tracklets, myAxisToSplit);
        }
        else {
            // use average like C linkTracklets
            pivot = (myUBounds[myAxisToSplit] + myLBounds[myAxisToSplit]) / 2.0;
        }

        // try to partition data
        for (uint i = 0; i < tracklets.size(); i++) {
            double val = tracklets[i].getPoint()[myAxisToSplit];
            
            if (val < pivot) {
                leftPointsAndValues.push_back(tracklets[i]);
            }
            else {
                rightPointsAndValues.push_back(tracklets[i]);
            }
        }
        
        // like in C linkTracklets, partition up data and if it
        // doesn't work well just partition arbitrarily...
        if ((leftPointsAndValues.size() == 0) || 
            (rightPointsAndValues.size() == 0)) {
            leftPointsAndValues.clear();
            rightPointsAndValues.clear();

            for (uint i = 0; i < tracklets.size(); i++) {
                if (i % 2 == 0) {
                    leftPointsAndValues.push_back(tracklets[i]);
                }
                else {
                    rightPointsAndValues.push_back(tracklets[i]);
                }
            }
        }

        nextAxis = (myAxisToSplit + 1) % (myK);
        
        TrackletTreeNode leftChild(allDetections,
                                   leftPointsAndValues,
                                   allLinkages,
                                   positionalErrorRa, positionalErrorDec,
                                   maxLeafSize, nextAxis, widths, lastId, 
                                   useMedian, splitWidest);
        
        myChildren.push_back(leftChild);
        
        TrackletTreeNode rightChild(allDetections,
                                    rightPointsAndValues, 
                                    allLinkages,
                                    positionalErrorRa, positionalErrorDec,
                                    maxLeafSize, nextAxis, widths, lastId, 
                                    useMedian, splitWidest);

        myChildren.push_back(rightChild);

        // get objects from child nodes
        for (unsigned int i = 0; i < myChildren.size(); i++) {
            const std::set<int> *childObjs = 
                (myChildren.at(i).getMyObjects());
            std::set<int>::const_iterator o;
            for (o = childObjs->begin(); o != childObjs->end(); o++) {
                myObjects.insert(*o);
            }
        }
    }


    // account for error in child tracklets.
    recalculateBoundsWithError(positionalErrorRa, positionalErrorDec);

}


const std::set<int> * TrackletTreeNode::getMyObjects()
{
    // this is populated at construction-time
    return & myObjects;
}


        

// these are to be used by linkTracklets.
const unsigned int TrackletTreeNode::getNumVisits() const
{
    return numVisits;
}


void TrackletTreeNode::addVisit() 
{
    numVisits++;
}


bool TrackletTreeNode::hasTracklet(unsigned int t)
{
    if (isLeaf()) {
        for (uint i = 0; i < myData.size(); i++) {
            if (myData[i].getValue() == t) {
                return true;
            }
        }
        return false;
    }
    else {
        if (hasLeftChild()) {
            if (getLeftChild()->hasTracklet(t)) {
                return true;
            }
        }
        if (hasRightChild()) {
            if (getRightChild()->hasTracklet(t)) {
                return true;
            }
        }
        return false;
    }
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
    






    



     
void TrackletTreeNode::recalculateBoundsWithError(
    double positionalErrorRa,
    double positionalErrorDec)
{
    if (isLeaf()) {
        
        // extend UBounds, LBounds with knowledge of velocity error.
        if (myData.size() <= 0) {
            LSST_EXCEPT(ProgrammerErrorException, 
      "Unexpected condition: tree node is leaf, but has no data.");
        }

        if (myK == 4) {
	
            // extend UBounds, LBounds in RA, Dec by position error.
            myUBounds.at(0) += positionalErrorRa;
            myLBounds.at(0) -= positionalErrorRa;
            myUBounds.at(1) += positionalErrorDec;
            myLBounds.at(1) -= positionalErrorDec;
            
            // find min/max RA, Dec velocities after accounting for error.
            
            for (unsigned int i = 0; i < myData.size(); i++) {
                std::vector<double> trackletPoint = 
                    myData.at(i).getPoint();
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
            /* this isn't really a tracklet, it's a track! */
            // extend UBounds, LBounds in RA, Dec by position error.
            myUBounds.at(0) += positionalErrorRa;
            myLBounds.at(0) -= positionalErrorRa;
            myUBounds.at(1) += positionalErrorDec;
            myLBounds.at(1) -= positionalErrorDec;
            // TBD: Hopefully Tim will write something smarter here,
            // for now used dummy values more or less.
            double DUMMY_ACCEL_ERROR = .01;
            double DUMMY_VEL_ERROR = .01;

            myUBounds.at(2) += DUMMY_VEL_ERROR;
            myLBounds.at(2) -= DUMMY_VEL_ERROR;
            myUBounds.at(3) += DUMMY_VEL_ERROR;
            myLBounds.at(3) -= DUMMY_VEL_ERROR;


            myUBounds.at(4) += DUMMY_ACCEL_ERROR;
            myLBounds.at(4) -= DUMMY_ACCEL_ERROR;
            myUBounds.at(5) += DUMMY_ACCEL_ERROR;
            myLBounds.at(5) -= DUMMY_ACCEL_ERROR;

            // guess that velocity error is basically just like in the
            // vel-only case. probably not. TBD, don't do something
            // this stupid.

            for (unsigned int i = 0; i < myData.size(); i++) {
                std::vector<double> trackletPoint = 
                    myData.at(i).getPoint();
                double thisDt = trackletPoint.at(6);

                double trackletRaV  = trackletPoint.at(2);
                double trackletDecV = trackletPoint.at(3);

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
    }
    else {
        // traverse children, then extend your own bounds.
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

        

}}//close lsst::Mops

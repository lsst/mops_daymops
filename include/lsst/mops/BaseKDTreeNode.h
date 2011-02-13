// -*- LSST-C++ -*-


#ifndef BASE_LSST_KDTREE_NODE_H
#define BASE_LSST_KDTREE_NODE_H

#include <iostream>

#include "lsst/mops/common.h"
#include "lsst/mops/PointAndValue.h"
#include "Exceptions.h"


/*
 * jmyers feb 2011
 * 
 * We really want two classes: KDTree and TrackletTree. KDTree supports
 * searching and has non-overlapping nodes and TrackletTree does not support
 * searching but does allow overlapping nodes.
 *
 * BaseKDTree will implement most of the stuff we need for both classes and they
 * will derive from it. Don't use this class on its own!
 * 
 */

namespace lsst {
namespace mops {

    
    template <class T, class RecursiveT>
    class BaseKDTreeNode {
    public: 
        /* create and populate a KDTreeNode.  Caller is trusted that the
         * Ubounds, LBounds provided correspond to the upper and lower
         * bounds of data in each dimension from pointsAndValues.
         *
         * automatically sets ref count to 1.
         */
        BaseKDTreeNode(
            const std::vector<PointAndValue <T> > &pointsAndValues, 
            unsigned int k, unsigned int maxLeafSize, 
            unsigned int myAxisToSplit, 
            const std::vector<double> &Ubounds,
            const std::vector<double> &LBounds, 
            unsigned int &lastId=0);

        /*
         * this addReference operator is to be used by KDTree *ONLY*.
         * This way multiple KDTrees can share the same nodes without
         * expensive copying (and since trees are static, this is fine.)
         *
         * it is up to KDTree to handle the creation and deletion of
         * these nodes based on refcount.
         */
        void addReference();

        void removeReference();
        
        unsigned int getRefCount();
        
        void debugPrint(int depth) const;
    
        const unsigned int getId() const;

        /* 
         * helper functions, should probably not be public... TBD.
         *
         *for a collection of points and values, find the median value
         * along the given axis and return it
         *
         * ASSUMES all points have size > axis 
        */
        double getMedianByAxis(std::vector<PointAndValue<T> > pointsAndValues,
                               unsigned int axis);

        double maxByAxis(std::vector<PointAndValue<T> > pointsAndValues,
                         unsigned int axis);


        const std::vector<double> *getUBounds() const;
        const std::vector<double> *getLBounds() const;

        
    protected:
        std::vector <RecursiveT > myChildren;        
        unsigned int myRefCount;
        unsigned int myK;
        std::vector <double> myUBounds;
        std::vector <double> myLBounds;
        std::vector <PointAndValue <T> > myData;
        unsigned int id;
    };




template <class T, class RecursiveT>
const std::vector<double> *BaseKDTreeNode<T, RecursiveT>::getUBounds() 
    const
{
    return &myUBounds;
}


template <class T, class RecursiveT>
const std::vector<double> *BaseKDTreeNode<T, RecursiveT>::getLBounds() 
    const
{
    return &myLBounds;
}


template <class T, class RecursiveT>
void BaseKDTreeNode<T, RecursiveT>::addReference()
{
    myRefCount++;
}




template <class T, class RecursiveT>
void BaseKDTreeNode<T, RecursiveT>::removeReference()
{
    myRefCount--;
}



template <class T, class RecursiveT>
const unsigned int BaseKDTreeNode<T, RecursiveT >::getId() const
{
    return id;
}



template <class T, class RecursiveT>
unsigned int BaseKDTreeNode<T, RecursiveT>::getRefCount()
{
    return myRefCount;
}





template <class T, class RecursiveT>
BaseKDTreeNode<T, RecursiveT>::BaseKDTreeNode(
    const std::vector<PointAndValue <T> > &pointsAndValues, 
    unsigned int k, 
    unsigned int maxLeafSize,       
    unsigned int myAxisToSplit, 
    const std::vector<double> & UBounds, 
    const std::vector<double> & LBounds,
    unsigned int &lastId)
{

    myRefCount = 1;
    myK = k;  
    myUBounds = UBounds;
    myLBounds = LBounds;
    lastId++;
    id = lastId;


    std::vector<double> rightChildUBounds, rightChildLBounds,     
        leftChildUBounds, leftChildLBounds;
    std::vector<PointAndValue <T> > leftPointsAndValues, rightPointsAndValues;
    double tmpMedian;

    typename std::vector<PointAndValue<T>,
        std::allocator<PointAndValue<T> > >::iterator myIter;

    PointAndValue<T> tmpPointAndValue;
    unsigned int nextAxis;

    if (myAxisToSplit >= myK){
        throw LSST_EXCEPT(BadParameterException, 
                          "BaseKDTreeNode: failed sanity check: asked to split data along dimension greater \
than dimensions of data\n");
    }

    if (myUBounds.size() != myLBounds.size()){
        throw LSST_EXCEPT(ProgrammerErrorException, 
                          "BaseKDTreeNode: failed sanity check: UBounds/LBounds vectors have different size\n");
    }


    /* 
       if myUBounds[i] == myLBounds[i] for i == 0,..k, but
       pointsAndValues.size() > 1, then ALL POINTS ARE EQUAL.
     
       this means that we need to be a leaf no matter what.
     
       for real-world (LSST) purposes, this should virtually NEVER
       happen, but it is theoretically possible, especially if
       this code is used for discrete values.
    */
    bool forceLeaf = true;
    for (unsigned int i = 0; 
         (i < myUBounds.size()) && (forceLeaf == true);
         i++) 
    {
        if (!areEqual(myUBounds[i], myLBounds[i])) {
            forceLeaf = false;
        } 
    }


    tmpMedian = getMedianByAxis(pointsAndValues, myAxisToSplit);
    /* 
     * catch the special case John Dailey from WISE found: if the
     * median value is also the max value, we can't just give all
     * values <= the median to the left child; this would be *all*
     * the values and can lead to infinite recursion.
     */
    if (tmpMedian == maxByAxis(pointsAndValues, myAxisToSplit)) {
        forceLeaf = true;
        /* in the future, we may consider changing this behavior;
         * rather than force a leaf node at this location, use <
         * rather than <= in choosing values to send to the left
         * child. */
    }



    /* case 1: this is a leaf. */
    if ((forceLeaf == true) || (pointsAndValues.size() <= maxLeafSize)) {
        myData = pointsAndValues;
    }
    else {
        /* case 2: not a leaf, so split data and into two chunks,
         * recursively create children with them by partitioning the data
         * along myAxisToSplit.  */


        rightChildUBounds = UBounds;

        rightChildLBounds = LBounds;
        rightChildLBounds[myAxisToSplit] = tmpMedian;

        leftChildUBounds = UBounds;
        leftChildUBounds[myAxisToSplit] = tmpMedian;

        leftChildLBounds = LBounds;
    
        //partition data for children
        for (unsigned int i = 0; i < pointsAndValues.size(); i++) {

            tmpPointAndValue = pointsAndValues[i];

            if (tmpPointAndValue.getPoint()[myAxisToSplit] <= tmpMedian) { 
                leftPointsAndValues.push_back(tmpPointAndValue);
            }
            else {
                rightPointsAndValues.push_back(tmpPointAndValue);
            }
      
            /* it would be really nice to delete
             * pointsAndValues[i] at this point, but since it
             * holds statically allocated data that's not really
             * possible... this calls for a possible redesign
             * (ugh)
             */
        }
    
        nextAxis = (myAxisToSplit + 1) % (myK);
    
        RecursiveT leftChild(leftPointsAndValues, k, maxLeafSize,
                             nextAxis, 
                             leftChildUBounds, leftChildLBounds, lastId);
    
        myChildren.push_back(leftChild);
    
        RecursiveT rightChild(rightPointsAndValues, k, maxLeafSize,
                              nextAxis, rightChildUBounds, 
                              rightChildLBounds, lastId);

        myChildren.push_back(rightChild);
    }
}




template <class T, class RecursiveT>
void BaseKDTreeNode<T, RecursiveT>::debugPrint(int depth) const
{
    std::cout << "KDTREENODE: Depth "<< depth << std::endl;;
    std::cout << "\tmy dims is "<< myK << std::endl;
    std::vector<double>::iterator myIter;
    
    std::cout << "UBounds: ";
    printDoubleVec(myUBounds);
    
    std::cout << std::endl << "LBounds: ";
    printDoubleVec(myLBounds);
    
    if (myChildren.size() == 0) { 
        
        // we're a leaf node
        typename std::vector<PointAndValue<T>,std::allocator<PointAndValue<T> > >::iterator dataIter;
        for (dataIter = myData.begin(); dataIter != myData.end(); dataIter++)
        {
            std::cout << "Data point: ";
            printDoubleVec(dataIter->getPoint());
        }
        std::cout << std::endl;
    }
    else {
        //we're a non-leaf, print the children's contents
        
        myChildren[0].debugPrint(depth + 1);
        myChildren[1].debugPrint(depth + 1);    
    }
}
    
    
template <class T, class RecursiveT>
double BaseKDTreeNode<T, RecursiveT>::maxByAxis(
    std::vector<PointAndValue<T> > pointsAndValues, 
    unsigned int axis) 
{
    double tmpMax = 0;
    bool foundMax = false;
    for (unsigned int i = 0; i < pointsAndValues.size(); i++) {
        double val = pointsAndValues.at(i).getPoint().at(axis);
        if ((val > tmpMax) || (foundMax == false)) {
            tmpMax = val;
            foundMax = true;
        }
    }
    return tmpMax;
}


 
template <class T, class RecursiveT>
double BaseKDTreeNode<T, RecursiveT>::getMedianByAxis(
    std::vector<PointAndValue<T> > pointsAndValues,
    unsigned int axis)
{
    double tmpMedian;
    std::vector<double> splitAxisPointData;
    typename std::vector<PointAndValue<T>,
        std::allocator<PointAndValue<T> > >::iterator myIter;
    
    for (myIter = pointsAndValues.begin(); 
         myIter != pointsAndValues.end();
         myIter++) {
        
        splitAxisPointData.push_back( myIter->getPoint()[axis]);
        
    }
    
    /*we know that we have at least one element here, so we're fine
     * (fastMedian demands non-empty vector)*/
    tmpMedian = fastMedian(splitAxisPointData);
  
    return tmpMedian;
    
}






}} // close namespace lsst::mops

#endif

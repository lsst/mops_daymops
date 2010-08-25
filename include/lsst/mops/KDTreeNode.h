// -*- LSST-C++ -*-


#ifndef LSST_KDTREE_NODE_H
#define LSST_KDTREE_NODE_H

#include <iostream>

#include "Exceptions.h"


/*
 * the actual implementation of nodes in the KDTree.
 * 
 *
 * Each node represents a partition of the K-dimensional space in
 * which our points lie.  If it is a leaf, it holds a vector of
 * points; if it is not a leaf, it holds a vector of children.  The
 * children represent a partitioning of the space along axis A s.t.
 * the left child (the 0 element) represents all data with value < a
 * along axis A, the right child (the 1 element) represents all data
 * with value >= a in axis A.
 *
 * In turn the children, if non-leaf nodes, will have their children
 * represent partitions of the space along axis B, where B is either
 * dimension A + 1 or dimension 1 if A == K.
 * 
 * THIS CLASS IS ONLY TO BE USED BY KDTree, and should *NEVER* be used elsewhere.
 * 
 */

namespace lsst {
namespace mops {


    template <class T>
    class KDTreeNode {
    public: 
        /* create and populate a KDTreeNode.  Caller is trusted that the
         * Ubounds, LBounds provided correspond to the upper and lower
         * bounds of data in each dimension from pointsAndValues.
         *
         * automatically sets ref count to 1.
         */
        KDTreeNode(std::vector<PointAndValue <T> > pointsAndValues, 
                   unsigned int k, unsigned int maxLeafSize, 
                   unsigned int myAxisToSplit, std::vector<double> Ubounds,
                   std::vector<double>LBounds, unsigned int &lastId=0);

        /*
         * this addReference operator is to be used by KDTree *ONLY*.
         * This way multiple KDTrees can share the same nodes without
         * expensive copying (and since trees are static, this is fine.)
         *
         * it is up to KDTree to handle the creation and deletion of these nodes
         * based on refcount.
         */
        void addReference();

        void removeReference();

        unsigned int getRefCount();

        std::vector<PointAndValue <T> > rangeSearch(std::vector<double> queryPt, 
                                                    double queryRange) const; 
    
        std::vector<PointAndValue <T> > 
        hyperRectangleSearch(const std::vector<double> &queryPt, 
                             const std::vector<double> &tolerances, 
                             const std::vector<GeometryType> &spaceTypesByDimension) const;

        // these are to be used by linkTracklets.
        bool hasLeftChild() const;
        bool hasRightChild() const;
        KDTreeNode<T> * getLeftChild();
        KDTreeNode<T> * getRightChild();
        /*return a pointer to a const vector of the per-axis upper bounds of
          this tree node.  The format is identical to the double vector in the
          pointsAndValues used to create this node.
         */
        const std::vector<double> *getUBounds() const;
        // same, but with the lower bounds rather than upper bounds.
        const std::vector<double> *getLBounds() const;
        const std::vector<PointAndValue <T> > * getMyData() const;
        bool isLeaf() const;
        
        void debugPrint(int depth) const;
    
        std::vector <KDTreeNode<T> > myChildren;
        
        const unsigned int getNumVisits() const;
        const unsigned int getId() const;
        void addVisit();
        
    private:
        /* for a collection of points and values, find the median value
           along the given axis and return it

           ASSUMES all points have size > axis 
        */
        double getMedianByAxis(std::vector<PointAndValue<T> > pointsAndValues,
                               unsigned int axis);

        double maxByAxis(std::vector<PointAndValue<T> > pointsAndValues,
                         unsigned int axis);
        unsigned int myRefCount;
        unsigned int myK;
        std::vector <double> myUBounds;
        std::vector <double> myLBounds;
        std::vector <PointAndValue <T> > myData;
        //these are used by linkTracklets
        unsigned int numVisits;
        unsigned int id;

    };



    // these are to be used by linkTracklets.
    template <class T>
    const unsigned int KDTreeNode<T>::getNumVisits() const
    {
        return numVisits;
    }

    template <class T>
    const unsigned int KDTreeNode<T>::getId() const
    {
        return id;
    }

    template <class T>
    void KDTreeNode<T>::addVisit() 
    {
        numVisits++;
    }


    template <class T>
    bool KDTreeNode<T>::hasLeftChild() const
    {
        if (myChildren.size() >= 1) {
            return true;
        }
        return false;
    }


    template <class T>
    bool KDTreeNode<T>::hasRightChild() const
    {
        if (myChildren.size() >= 2) {
            return true;
        }
        return false;
    }


    template <class T>
    KDTreeNode<T> * KDTreeNode<T>::getLeftChild()
    {
        if (!hasLeftChild()) {
            return NULL;
        }
        return &(myChildren.at(0));
    }

 
    template <class T>
    KDTreeNode<T> * KDTreeNode<T>::getRightChild()
    {
        if (!hasRightChild()) {
            return NULL;
        }
        return &(myChildren.at(1));
    }


    template <class T>
    const std::vector<double> *KDTreeNode<T>::getUBounds() const
    {
        return &myUBounds;
    }


    template <class T>
    const std::vector<double> *KDTreeNode<T>::getLBounds() const
    {
        return &myLBounds;
    }


    template <class T>
    const std::vector<PointAndValue <T> > * KDTreeNode<T>::getMyData() const
    {
        if (!isLeaf()) {
            LSST_EXCEPT(BadParameterException,
                        "KDTreeNode got request for data, but is not a leaf node.");
        }
        return &myData;
    }


    template <class T>
    bool KDTreeNode<T>::isLeaf() const
    {
        //std::cout << "KDTreeNode.h: myChildren.size() == " << myChildren.size() << std::endl;
        if (myChildren.size() == 0) 
        {
            return true;
        }
        return false;
    }
    




    template <class T>
    void KDTreeNode<T>::addReference()
    {
        myRefCount++;
    }

    template <class T>
    void KDTreeNode<T>::removeReference()
    {
        myRefCount--;
    }

    template <class T>
    unsigned int KDTreeNode<T>::getRefCount()
    {
        return myRefCount;
    }



    template <class T>
    KDTreeNode<T>::KDTreeNode(std::vector<PointAndValue <T> > pointsAndValues, 
                              unsigned int k, unsigned int maxLeafSize,       
                              unsigned int myAxisToSplit, std::vector<double> UBounds, 
                              std::vector<double> LBounds, unsigned int &lastId)
    {
        myRefCount = 1;
        myK = k;  
        myUBounds = UBounds;
        myLBounds = LBounds;
       
        lastId++;
        id = lastId;

        numVisits = 0;

        std::vector<double> rightChildUBounds, rightChildLBounds,     
            leftChildUBounds, leftChildLBounds;
        std::vector<PointAndValue <T> > leftPointsAndValues, rightPointsAndValues;
        double tmpMedian;
        typename std::vector<PointAndValue<T>,std::allocator<PointAndValue<T> > >::iterator myIter;
        PointAndValue<T> tmpPointAndValue;
        unsigned int nextAxis;

        if (myAxisToSplit >= myK){
            throw LSST_EXCEPT(BadParameterException, 
                              "KDTreeNode: failed sanity check: asked to split data along dimension greater than dimensions of data\n");
        }

        if (myUBounds.size() != myLBounds.size()){
            throw LSST_EXCEPT(ProgrammerErrorException, 
                              "KDTreeNode: failed sanity check: UBounds/LBounds vectors have different size\n");
        }


        /* 
           if myUBounds[i] == myLBounds[i] for i == 0,..k, but pointsAndValues.size() > 1, then
           ALL POINTS ARE EQUAL.
     
           this means that we need to be a leaf no matter what.
     
           for real-world (LSST) purposes, this should virtually NEVER happen, but it is theoretically
           possible, especially if this code is used for discrete values.
        */
        bool forceLeaf = true;
        for (unsigned int i = 0; (i < myUBounds.size()) && (forceLeaf == true); i++) 
        {
            if (!areEqual(myUBounds[i], myLBounds[i])) {
                forceLeaf = false;
            } 
        }


        tmpMedian = getMedianByAxis(pointsAndValues, myAxisToSplit);
        /* 
         * catch the special case John Dailey from WISE found: if the median value is
         * also the max value, we can't just give all values <= the median to the left
         * child; this would be *all* the values and can lead to infinite recursion.
         */
        if (tmpMedian == maxByAxis(pointsAndValues, myAxisToSplit)) {
            forceLeaf = true;
            /* in the future, we may consider changing this behavior; rather than
             * force a leaf node at this location, use < rather than <= in choosing
             * values to send to the left child. */
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
      
                /* it would be really nice to delete pointsAndValues[i] at this
                 * point, but since it holds statically allocated data that's
                 * not really possible... this calls for a possible redesign (ugh)
                 */
            }
    
            nextAxis = (myAxisToSplit + 1) % (myK);
    
            KDTreeNode leftChild(leftPointsAndValues, k, maxLeafSize,
                                 nextAxis, leftChildUBounds, leftChildLBounds, lastId);
    
            myChildren.push_back(leftChild);
    
            KDTreeNode rightChild(rightPointsAndValues, k, maxLeafSize,
                                  nextAxis, rightChildUBounds, rightChildLBounds, lastId);

            myChildren.push_back(rightChild);
        }
    }




    template <class T>
    void KDTreeNode<T>::debugPrint(int depth) const
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


    template <class T>
    double KDTreeNode<T>::maxByAxis(std::vector<PointAndValue<T> > pointsAndValues, 
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


    template <class T>
    double KDTreeNode<T>::getMedianByAxis(std::vector<PointAndValue<T> > pointsAndValues,
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
  
    /*we know that we have at least one element here, so we're fine (fastMedian demands non-empty vector)*/
    tmpMedian = fastMedian(splitAxisPointData);
  
    return tmpMedian;

}






    template <class T>
    std::vector<PointAndValue <T> > KDTreeNode<T>::rangeSearch(std::vector<double> queryPt, \
                                                               double queryRange) const
    {
        /* if we are not within queryRange[i] of queryPt[i] on either edge,
           and in any direction, then we cannot possibly be within range of
           the search; just return.

           otherwise, pass the buck to the kids if a non-leaf or do an
           actual search of the data if we are a leaf.
        */

        std::vector<PointAndValue <T> > myResults;
        bool isInRange = true;
        std::vector<double> query(1);
        std::vector<double> uBound(1);
        std::vector<double> lBound(1);
        typename std::vector<PointAndValue<T>,std::allocator<PointAndValue<T> > >::iterator dataIter;

        for (unsigned int i = 0; i < myK; i++)
        {
            query[0] = queryPt[i];
            uBound[0] = myUBounds[i];
            lBound[0] = myLBounds[i];
            if ((euclideanDistance(query, uBound, 1) > queryRange) && 
                (euclideanDistance(query, lBound, 1) > queryRange))
                isInRange = false;
        }
        if (isInRange == true)
        {
            if (myChildren.size() == 0)
            {
                /* this is a leaf node, search through the data */
                for (dataIter = myData.begin(); dataIter != myData.end();
                     dataIter++) {
                    if (euclideanDistance(queryPt, dataIter->getPoint(), myK) <= queryRange)
                    {
                        myResults.push_back(*dataIter);
                    }
                }
            }
            else {
                /* not a leaf node, so just pass the buck */

                std::vector<PointAndValue<T> > childResults;
                for (unsigned int i = 0; i < 2; i++) {
                    childResults = myChildren[i].rangeSearch(queryPt, queryRange);
                    for (dataIter = childResults.begin(); dataIter != childResults.end();
                         dataIter++) {
                        myResults.push_back(*dataIter);
                    }
                }
            }
        }
        return myResults;
    }






template <class T>
std::vector<PointAndValue <T> > 
KDTreeNode<T>::hyperRectangleSearch(const std::vector<double> &queryPt,
				    const std::vector<double> &tolerances,
				    const std::vector<GeometryType> &spaceTypesByDimension) const 
{
    std::vector<PointAndValue <T> > myResults;

    /* 
     * just like for rangeSearch, we want to return nothing if queryPt is too
     * far from our representative space; otherwise, we want to either pass the
     * buck to our children and accumulate their results, or we want to search
     * the data ourselves.
     */

    bool isInRange = true;

    for (unsigned int i = 0; i < spaceTypesByDimension.size(); i++) {
        if (spaceTypesByDimension[i] == CIRCULAR_DEGREES) {
            if ((myUBounds[i] != convertToStandardDegrees(myUBounds[i]))
                ||
                (myLBounds[i] != convertToStandardDegrees(myLBounds[i]))
                ||
                (queryPt[i] != convertToStandardDegrees(queryPt[i]))) {
                std::cerr << "Requested values: " << myUBounds[i] << ", " 
                          << myLBounds[i] << " " << queryPt[i] << " and i == " << i << std::endl;
                throw LSST_EXCEPT(BadParameterException,  
                                  "KDTreeNode: Data error: got that dimension is of type CIRCULAR_DEGREES but data and/or query do not lie along [0,360).");
            }
        }
    }
    
    for (unsigned int i = 0; (i < myK) && (isInRange == true); i++) {
        double UBoundTolerance = queryPt[i] + tolerances[i];
        double LBoundTolerance = queryPt[i] - tolerances[i];
        /* NB: since tolerances may be up to 180 degrees, it is necessary to do
         *two* tests here - one for each semicircle, one for the
         'bottom' half.  this is due to the slightly curious implementation of regionsOverlap1D*/
        if ((regionsOverlap1D(myLBounds[i], myUBounds[i], 
                                      UBoundTolerance, queryPt[i],
                                      spaceTypesByDimension[i]) == false) && 
	    (regionsOverlap1D(myLBounds[i], myUBounds[i], 
				      queryPt[i], LBoundTolerance,
				      spaceTypesByDimension[i]) == false)) { 
            isInRange = false;
        }
    }

    if (isInRange == true) {
        if (myChildren.size() != 0) {
            /* punt to the children */
            for (unsigned int i = 0; i < myChildren.size(); i++) {
                std::vector<PointAndValue <T> > childResults;
                childResults = myChildren[i].hyperRectangleSearch(queryPt, tolerances, 
                                                                  spaceTypesByDimension);

                for (unsigned int j = 0; j < childResults.size(); j++) {                    

                    myResults.push_back(childResults[j]);
                }                   
            }

        }
        else { 
            /* do the actual searching */
            for (unsigned int i = 0; i < myData.size(); i++) {                
                isInRange = true;
                for (unsigned int j = 0; j < myK; j++) {
                    std::vector<double> point = myData[i].getPoint();
                    if (distance1D(point[j], queryPt[j], spaceTypesByDimension[j]) > tolerances[j]) {
                        isInRange = false;
                    }
                }
                if (isInRange == true) {
                    myResults.push_back(myData[i]);
                }
            }
            
        }
        
    }
    
    return myResults;
}

 





}} // close namespace lsst::mops

#endif

// -*- LSST-C++ -*-


#ifndef LSST_KDTREE_NODE_H
#define LSST_KDTREE_NODE_H

#include <iostream>

#include "lsst/mops/PointAndValue.h"
#include "lsst/mops/BaseKDTreeNode.h"
#include "lsst/mops/Exceptions.h"

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

static unsigned int lastId2_NULL = 0;

namespace lsst {
namespace mops {


    template <class T>
    class KDTreeNode: public BaseKDTreeNode <T, KDTreeNode<T> > {
    public: 

        /* create and populate a KDTreeNode.  Caller is trusted that the
         * Ubounds, LBounds provided correspond to the upper and lower
         * bounds of data in each dimension from pointsAndValues.
         *
         * automatically sets ref count to 1.
         */
        KDTreeNode(const std::vector<PointAndValue <T> > &pointsAndValues, 
                   unsigned int k, 
                   unsigned int maxLeafSize, 
                   unsigned int myAxisToSplit, 
                   const std::vector<double> &Ubounds,
                   const std::vector<double> &LBounds, 
                   unsigned int &lastId=lastId2_NULL) 
            
            : BaseKDTreeNode<T, KDTreeNode<T> >(pointsAndValues, 
                                                k, maxLeafSize, 
                                                myAxisToSplit, Ubounds,
                                                LBounds, lastId) {}
        


        std::vector<PointAndValue <T> > rangeSearch(
            std::vector<double> queryPt, 
            double queryRange) const; 
    
        std::vector<PointAndValue <T> > 
        hyperRectangleSearch(const std::vector<double> &queryPt, 
                             const std::vector<double> &tolerances, 
                             const std::vector<GeometryType> &spaceTypesByDimension) const;

        

    };







    template <class T>
    std::vector<PointAndValue <T> > 
    KDTreeNode<T>::rangeSearch(std::vector<double> queryPt,
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

        for (unsigned int i = 0; i < this->myK; i++)
        {
            query[0] = queryPt[i];
            uBound[0] = this->myUBounds[i];
            lBound[0] = this->myLBounds[i];
            if ((euclideanDistance(query, uBound, 1) > queryRange) && 
                (euclideanDistance(query, lBound, 1) > queryRange))
                isInRange = false;
        }
        if (isInRange == true)
        {
            if (this->myChildren.size() == 0)
            {
                /* this is a leaf node, search through the data */
                for (dataIter = this->myData.begin(); dataIter != this->myData.end();
                     dataIter++) {
                    if (euclideanDistance(queryPt, dataIter->getPoint(), this->myK) <= queryRange)
                    {
                        myResults.push_back(*dataIter);
                    }
                }
            }
            else {
                /* not a leaf node, so just pass the buck */

                std::vector<PointAndValue<T> > childResults;
                for (unsigned int i = 0; i < 2; i++) {
                    childResults = this->myChildren[i].rangeSearch(queryPt, queryRange);
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
				    const std::vector<GeometryType> &spaceTypesByDimension) 
    const 
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
            if ((this->myUBounds[i] != convertToStandardDegrees(this->myUBounds[i]))
                ||
                (this->myLBounds[i] != convertToStandardDegrees(this->myLBounds[i]))
                ||
                (queryPt[i] != convertToStandardDegrees(queryPt[i]))) {
                std::cerr << "Requested values: " << this->myUBounds[i] << ", " 
                          << this->myLBounds[i] << " " << queryPt[i] 
                          << " and i == " << i << std::endl;
                throw LSST_EXCEPT(BadParameterException,  
                                  "KDTreeNode: Data error: got that dimension is of type CIRCULAR_DEGREES but data and/or query do not lie along [0,360).");
            }
        }
    }
    
    for (unsigned int i = 0; (i < this->myK) && (isInRange == true); i++) {
        double UBoundTolerance = queryPt[i] + tolerances[i];
        double LBoundTolerance = queryPt[i] - tolerances[i];
        /* NB: since tolerances may be up to 180 degrees, it is necessary to do
         *two* tests here - one for each semicircle, one for the
         'bottom' half.  this is due to the slightly curious implementation of regionsOverlap1D*/
        if ((regionsOverlap1D(this->myLBounds[i], this->myUBounds[i], 
                              UBoundTolerance, queryPt[i],
                              spaceTypesByDimension[i]) == false) && 
	    (regionsOverlap1D(this->myLBounds[i], this->myUBounds[i], 
                              queryPt[i], LBoundTolerance,
                              spaceTypesByDimension[i]) == false)) { 
            isInRange = false;
        }
    }

    if (isInRange == true) {
        if (this->myChildren.size() != 0) {
            /* punt to the children */
            for (unsigned int i = 0; i < this->myChildren.size(); i++) {
                std::vector<PointAndValue <T> > childResults;
                childResults = this->myChildren[i].hyperRectangleSearch(queryPt, tolerances, 
                                                                  spaceTypesByDimension);

                for (unsigned int j = 0; j < childResults.size(); j++) {                    

                    myResults.push_back(childResults[j]);
                }                   
            }

        }
        else { 
            /* do the actual searching */
            for (unsigned int i = 0; i < this->myData.size(); i++) {                
                isInRange = true;
                for (unsigned int j = 0; j < this->myK; j++) {
                    std::vector<double> point = this->myData[i].getPoint();
                    if (distance1D(point[j], queryPt[j], spaceTypesByDimension[j]) > tolerances[j]) {
                        isInRange = false;
                    }
                }
                if (isInRange == true) {
                    myResults.push_back(this->myData[i]);
                }
            }
            
        }
        
    }
    
    return myResults;
}

 





}} // close namespace lsst::mops

#endif

// -*- LSST-C++ -*-



/* jmyers 7/25/08
 *
 * this file is the default interface to our KDTree implementation.  There is
 * also a KDTreeNode class, which implements most of the work, but KDTree holds
 * onto the root of the tree and serves as the interface for building,
 * searching, and destroying trees.
 * 
 * A KDTree is an K-dimensional spatial data structure which allows for fast
 * spatial searches, but has no balancing algorithm.  As a result, data can only
 * be added to the tree once.
 *
 * The data type used by the KDTree is PointAndValue, which is a general-purpose
 * class which holds an K-dimensional vector of doubles (the spatial data, a
 * point) and some other data (the value)
 *
 * Since PointAndValue is a template class, KDTree also becomes a template
 * class. Template classes must keep all code in the header files.  This is a
 * really ugly and unfortunate mess.  To keep things readable, declarations will
 * exist in KDTree.cc, but including KDTree.h will cause KDTree.cc to be
 * included as well.
 */


#ifndef LSST_KDTREE_H
#define LSST_KDTREE_H

#include <vector>
#include <cmath>
#include <algorithm>

#include "common.h"  
#include "Exceptions.h"
#include "PointAndValue.h"
#include "BaseKDTree.h"
#include "KDTreeNode.h"
  

namespace lsst {
namespace mops {

    
    template <class T>
    class KDTree : public BaseKDTree<T, KDTree<T> > {
    public:


        /*
         * instantiate a KDTree with data from vector PointsAndValues (the
         * points are used to distribute data throughout the tree, the
         * associated values are some useful data you provide).  k is the
         * number of dimensionality of the space searched by the tree.
         * maxLeafSize is the maximum number of items in any leaf of the
         * tree, used for very particular optimizations.
         *
         * note that each element of pointsAndValues contains a vector of
         * N-dimensional 'points'; all vectors must have at least 'dimensions'
         * points where 'dimensions' is 1 or greater.  If points are of length >
         * dimensions, only the first dimensions elements are used.
         *
         * maxLeafSize should be a positive integer.
         *
         */
        KDTree(std::vector<PointAndValue <T> > pointsAndValues, 
               unsigned int k, unsigned int maxLeafSize);


        // public access to the build a tree after setting it up.
        void buildFromData(std::vector<PointAndValue <T> > pointsAndValues, 
                           unsigned int k, unsigned int maxLeafSize) 
            { 
                BaseKDTree<T, KDTree<T> >(pointsAndValues,
                                          k, maxLeafSize); 
            }

        /* rangeSearch: given a point queryPt and a range queryRange, treat all
         * the data points in the tree as though they live on the same
         * homogenous Euclidean space (i.e. all dimensions have the ame units)
         * and return all points within queryRange of queryPt.
         * 
         * if queryPt is not of the same dimensions as the tree, an exception will be thrown.
         */
        std::vector<PointAndValue <T> > 
        rangeSearch(std::vector<double> queryPt,
                    double queryRange) const;


        /*
         * RADecRangeSearch : Perform a range search on the RA and Dec axes of
         * some data, and optionally perform a hyperRectangle search on the
         * *other* axes (if any).
         *
         * The range search searches a range along the surface of the sphere
         * (i.e, the great-circle distance, not a euclidean approximation).  The
         * other dimensions are treated independently as either euclidean or
         * circular spaces.
         *
         * This, admittedly, has some horribly complex call semantics (sorry,
         * but making this code highly general required making the interface
         * complex.)
         *
         * This is probably best explained with examples.
         * 
         *
         * EXAMPLE 1:
         *
         * Basically, if you have a 2-D tree of RA and Dec specified in degrees,
         * and you would like to search 4 degrees around the RA, Dec point (180, -10),
         * call with the following parameters:
         * 
         * queryPt = [180, -10]
         * queryRange = 4.0
         * otherDimsPoint = []
         * otherTolerances = []
         * spaceTypesByDimension = [RA_DEGREES, DEC_DEGREES]
         * 
         * (note that otherTolerances is basically ignored).
         *
         *
         * EXAMPLE 2:
         * 
         * Say you are tracking moving objects which are moving on a spehere.
         * You may Have 4D data, in the form RA, Dec, velocity, and movement
         * angle.  
         * 
         * Say you'd like to find all points which are within 4.0 degrees of the
         * RA, Dec point [90, 70] on the sphere, *and* you would also like to
         * return only points which are moving at velocity 10 +/- 3 and at
         * the angle 80 degrees +/- 5 degrees.  Your parameters would be:
         * 
         * queryPt = [90,70]
         * queryRange = 4.0
         * otherDimsPoint = [ 10, 80 ]
         * otherTolerances = [ 3, 5 ]
         * spaceTypeByDimension = [RA_DEGREES, DEC_DEGREES, 
         *                         EUCLIDEAN, CIRCULAR_DEGREES]
         *
         * 
         * 
         * EXAMPLE 3:
         * 
         * Say your data was originally of the form [time in days, Dec, RA].
         * (Why you did it this way, with Dec first and RA second, I don't
         * know. That's your business.)  You would like to get points within 2
         * degrees around the RA,Dec point [180, -80] and within 3 days of day
         * 50.
         * 
         * queryPt = [180, -80]
         * queryRange = 2.0
         * otherDimsPoint = [50]
         * otherTolerances = [3]
         * spaceTypesByDimension = [ Euclidean, 
         *                           DEC_DEGREES, RA_DEGREES]
         *
         * 
         * More formal documentation is TBD.  Hopefully the examples are
         * sufficiently illustrative.
         */
        std::vector<PointAndValue <T> > 
            RADecRangeSearch(const std::vector<double> &RADecQueryPoint, 
                             double RADecQueryRange, 
			     const std::vector<double> &otherDimsPoint,
                             const std::vector<double> &otherDimsTolerances,
                             const std::vector<GeometryType> 
                                &spaceTypesByDimension) 
            const; 
    
        /* execute a hyperRectangle-shaped range search around queryPt.  queryPt and tolerances
         * must have the same dimensions as the tree.
         *
         * if centerPt is [ a, b, c, ... n] and 
         *
         * tolerances is [ x, y, z, .... m ]
         *
         *  any PointAndValue with point [q, r, s,... t] is returned iff
         *  
         *  a - x < q < a + x
         *  b - y < r < b + y
         *  c - z < s < c + z
         *  ...
         *  n - m < t < n + m
         *
         * spaceTypesByDimension should be a length-k vector describing the
         * types of spaces; for example, if searching data of the form
         * 
         * [x-position, y-position, direction angle in degrees, velocity]
         *
         * then you would want to use a vector like this:
         * { Euclidean, Euclidean, Circular_degrees, Euclidean }.
         *
         * However, if searching simply in [x-position, y-position, time ]
         * you would use { Euclidean, Euclidean, Euclidean. }
         *
         * Note that this style of searching is not yet supported by rangeSearch; 
         * currently all range searches (hypersphere searches) treat the space
         * as k-dimensional and Euclidean.
         *
         */
        std::vector<PointAndValue <T> > 
        hyperRectangleSearch(const std::vector<double> &queryPt, 
                             const std::vector<double> &tolerances, 
                             const std::vector<GeometryType> 
                                    &spaceTypesByDimension) 
            const;


    };









template <class T>
KDTree<T>::KDTree(std::vector<PointAndValue <T> > pointsAndValues,
                       unsigned int k, unsigned int maxLeafSize) 
{
    this->setUpEmptyTree();
    buildFromData(pointsAndValues, k, maxLeafSize);
}











template <class T>
std::vector<PointAndValue <T> > 
KDTree<T>::rangeSearch(std::vector<double> queryPt,
		       double queryRange) const
{
    /* sanity check */
    if (queryPt.size() != this->myK)
    {
        throw LSST_EXCEPT(BadParameterException, "KDTree::rangeSearch:  got myK != queryPoint size");
    }
    /* just punt to the KDTreeNode. */
    return this->myRoot->rangeSearch(queryPt, queryRange);
}
    
    
    
    
    
    

template <class T>
std::vector<PointAndValue <T> > 
KDTree<T>::RADecRangeSearch(const std::vector<double> &RADecQueryPoint, 
                            double RADecQueryRange, 
			    const std::vector<double> &otherDimsPoint,
                            const std::vector<double> &otherDimsTolerances,
                            const std::vector<GeometryType> &spaceTypesByDimension)  const
{
    /*
     * this function is implemented by finding a series of rectangles which will
     * enscribe the actual circle along the surface of the sphere. These
     * rectangles are searched with hyperRectangleSearch (along with the
     * additional parameters) and results are pruned.
     */

    // step 1: check input data.
    if ((RADecQueryPoint.size() != 2) || 
        (otherDimsPoint.size() != this->myK - 2) || 
        (otherDimsTolerances.size() != this->myK - 2) || 
        (spaceTypesByDimension.size() != this->myK) ||
        (RADecQueryRange <= 0.0))
    {
        throw LSST_EXCEPT(BadParameterException, 
    "KDTree::RADecRangeSearch called with illegal parameters.");
    }
    int RADimIndex = -1;
    int DecDimIndex = -1;
    for (unsigned int i = 0; i < spaceTypesByDimension.size(); i++) {
        if (spaceTypesByDimension.at(i) == RA_DEGREES) {
            RADimIndex = i;
        }
        else if (spaceTypesByDimension.at(i) == DEC_DEGREES) {
            DecDimIndex = i;
        }
    }
    if ((RADimIndex == -1) || (DecDimIndex == -1)) {
        LSST_EXCEPT(BadParameterException,
                    "KDTree::RADecRangeSearch called with spaceTypesByDimension missing either RA, Dec, or both - this is illegal");        
    }
    double RACenter =  convertToStandardDegrees(RADecQueryPoint.at(0));
    double DecCenter = convertToStandardDegrees(RADecQueryPoint.at(1));

    /* now, find a set of rectangles which enscribe the RA Dec range.
     * there will be at most 2.

     * simple case: the range does not pass over the north or south pole. You
     * get one rectangle, which enscribes the circle.
     * 
     * complicated case 1: the range passes over either the north or south pole.
     * you get one rectangle, which has RA width 360 (recall that at the pole,
     * 360 degrees in RA is basically an infinitely small area).
     *
     * Complicated case 2: The dec range passes over BOTH poles, meaning you
     * actually now have to do a brute force search over all the RA, Dec data! 
     */

    std::vector<double> realQueryPoint;
    std::vector<double> realQueryTolerances;
    std::vector<GeometryType> realQueryTypes;
    //Constants
    const double northPole_Dec = 90;
    const double southPole_Dec = 270;

    double RAHalfWidth = 0;
    double DecHalfWidth = 0;

    if ((circularShortestPathLen_Deg(DecCenter, northPole_Dec) < RADecQueryRange) ||
        (circularShortestPathLen_Deg(DecCenter, southPole_Dec) < RADecQueryRange))
    {
        // this query range crosses a pole, ergo a complicated case
        RAHalfWidth = 180; 
        /* at a pole, we need to search all 360 degrees around
         * the pole.
         */
        if ((circularShortestPathLen_Deg(DecCenter, northPole_Dec) < RADecQueryRange)  && 
            (circularShortestPathLen_Deg(DecCenter, southPole_Dec) < RADecQueryRange))
        {
            /* the query range crosses both poles - so we really search the whole sphere! */
            DecHalfWidth = 180;            
        }
        else {
            DecHalfWidth = RADecQueryRange;
        }
    }
    else { 
        // case 1: just one query 
        RAHalfWidth = maxOfTwo(
            arcToRA(DecCenter + RADecQueryRange, RADecQueryRange),
            arcToRA(DecCenter - RADecQueryRange, RADecQueryRange)
            );        
        DecHalfWidth = RADecQueryRange;
    }
    //build the real search params, pass them off to hyperRectangleSearch.
    unsigned int otherParamsIndexCounter = 0;
    for (unsigned int i = 0; i < this->myK; i++) {
        if (spaceTypesByDimension.at(i) == RA_DEGREES) {
            realQueryPoint.push_back(RACenter);
            realQueryTolerances.push_back(RAHalfWidth);            
            realQueryTypes.push_back(CIRCULAR_DEGREES);
        }
        else if (spaceTypesByDimension.at(i) == DEC_DEGREES) {
            realQueryPoint.push_back(DecCenter);
            realQueryTolerances.push_back(DecHalfWidth);
            realQueryTypes.push_back(CIRCULAR_DEGREES);
        }
        else {
            realQueryPoint.push_back(otherDimsPoint.at(otherParamsIndexCounter));
            realQueryTolerances.push_back(otherDimsTolerances.at(otherParamsIndexCounter));
            realQueryTypes.push_back(spaceTypesByDimension.at(i));
            otherParamsIndexCounter++;
        }
    }
             
    // now do a hyperRectangleSearch with this data. 
    std::vector<PointAndValue <T> > searchResults;
    searchResults = this->myRoot->hyperRectangleSearch(
        realQueryPoint, realQueryTolerances, realQueryTypes);    

    //prune results on angular distance around the center of the RA, Dec query.
    std::vector<PointAndValue <T> > prunedResults;
    for (unsigned int i = 0; i < searchResults.size(); i++) {
        std::vector<double> point = searchResults.at(i).getPoint(); 
        double resultRA = point.at(RADimIndex);
        double resultDec = point.at(DecDimIndex);

        if (angularDistanceRADec_deg(resultRA, resultDec, RACenter, DecCenter) 
            < RADecQueryRange) {
            prunedResults.push_back(searchResults.at(i));
        }
    }

    return prunedResults;
}
    



template <class T>
std::vector<PointAndValue <T> > 
KDTree<T>::hyperRectangleSearch(const std::vector<double> &queryPt,
				const std::vector<double> &tolerances,
				const std::vector<GeometryType> &spaceTypesByDimensions) const
{
    
  /* sanity check */
  if ((queryPt.size() != this->myK) || (tolerances.size() != this->myK) || 
      (spaceTypesByDimensions.size() != this->myK)) {
      throw LSST_EXCEPT(BadParameterException, 
 "EE: QueryPt must have dimensions at least equal to dimensions of tree.\n");
  }
  if (this->hasData != true) {
      // if we are queried, but do not have any data, return nothing.
      std::vector<PointAndValue <T> > toRet;
      toRet.clear();
      return toRet;
    }
  else {
      for (unsigned int i = 0; i < this->myK; i++) {
          
          if (((spaceTypesByDimensions[i] == CIRCULAR_DEGREES) 
               && (tolerances[i] > 180.)) ||
              (((spaceTypesByDimensions[i] == CIRCULAR_RADIANS) 
                && (tolerances[i] > M_PI)))) {
              throw LSST_EXCEPT(BadParameterException,
                                "EE: KDTree.hyperRectangleSearch: searching a radius greater than 180 degrees (pi radians) is meaningless.\n");
          }
      }
      
      /* just punt to the KDTreeNode. */
      
      return this->myRoot->hyperRectangleSearch(
          queryPt, tolerances, spaceTypesByDimensions);
  }
}







}} // close namespace lsst::mops

#endif

// -*- LSST-C++ -*-



/* jmyers 2/12/2011
 *
 * Where we once had the monolithic KDTree, there are now
 * 
 * BaseKDTree -> (derived) -> KDTree (for range searching)
 *    |
 *    \______-> (derived) -> TrackletTree (for manually walking 
 *                                          in linkTracklets.)
 * 
 * TreeNodes now follow a similar hierarchy. 
 *
 * This is the Base class, which manages memory issues (allocation/deletion,e
 * tc.) and functionality shared by both classes.
 */


#ifndef LSST_BASEKDTREE_H
#define LSST_BASEKDTREE_H

#include <vector>
#include <cmath>
#include <algorithm>

#include "common.h"  
#include "Exceptions.h"
#include "PointAndValue.h"
#include "BaseKDTreeNode.h"
  

namespace lsst {
namespace mops {

    
    template <class T, class TreeNodeClass>
    class BaseKDTree {
    public:

        /* default constructor - don't instantiate with any data.         
         */
        BaseKDTree();

        /* copy constructor */
        BaseKDTree(const BaseKDTree<T, TreeNodeClass> &source);
        
        void debugPrint() const;
        
        unsigned int size() const;

        BaseKDTree<T, TreeNodeClass>& operator=(
            const BaseKDTree<T, TreeNodeClass> &rhs);

        ~BaseKDTree();

        TreeNodeClass *myRoot;

    protected:

        /* 
         * populates the tree with given data.  See comments KDTree's second
         * constructor.
         */
        void buildFromData(std::vector<PointAndValue <T> > pointsAndValues, 
               unsigned int k, unsigned int maxLeafSize);

        void copyTree(const BaseKDTree<T, TreeNodeClass> &source);
        void setUpEmptyTree();
        void clearPrivateData();
        bool hasData;
        unsigned int myK;
        std::vector <double> myUBounds;
        std:: vector <double> myLBounds;

        unsigned int mySize;
    };







template <class T, class TreeNodeClass>
void BaseKDTree<T, TreeNodeClass>::clearPrivateData()
{
    if (hasData == true) {
        myRoot->removeReference();
        unsigned int refCount = myRoot->getRefCount();
        if (refCount == 0) {
            /*std::cerr << "DD: Deleting root " << std::endl;*/
            delete myRoot;
        }
        if (refCount < 0) {
            throw LSST_EXCEPT(ProgrammerErrorException,
              "EE: IMPOSSIBLE CASE: KDTreeNode has negative refcount");
        }
        hasData = false;
        myRoot = NULL;
    }
    myK = 0;
    myUBounds.clear();
    myLBounds.clear();
}
 




template <class T, class TreeNodeClass>
BaseKDTree<T, TreeNodeClass>::~BaseKDTree()
{
    /*std::cerr << "DD: Entering ~KDTree " << std::endl;*/
    clearPrivateData();
    /*std::cerr << "DD: exiting ~KDTree " << std::endl;*/
}





template <class T, class TreeNodeClass>
void BaseKDTree<T, TreeNodeClass>::setUpEmptyTree() 
{
    hasData = false;
    myRoot = NULL;
    myK = 0;    
    myUBounds.clear();
    myLBounds.clear();
}




template <class T, class TreeNodeClass>
BaseKDTree<T, TreeNodeClass>::BaseKDTree(
    const BaseKDTree<T, TreeNodeClass> &source) 
{
    setUpEmptyTree();
    /*std::cerr << "DD: entering copy constructor..." << std::endl;*/
    copyTree(source);
    /*std::cerr << "DD: exiting copy constructor..." << std::endl;*/
}



template<class T, class TreeNodeClass>
void BaseKDTree<T, TreeNodeClass>::copyTree(
    const BaseKDTree<T, TreeNodeClass> &source) 
{
    // check for self-assignment
    if (this != &source) { 
        clearPrivateData();
        myK = source.myK;
        if (source.hasData) {        
            myRoot = source.myRoot;
            myUBounds = source.myUBounds;
            myLBounds = source.myLBounds;
        /* make sure this node now knows it has an additional owner. */
            myRoot->addReference();
            hasData = true;
        }
    }
}


template <class T, class TreeNodeClass>
BaseKDTree<T, TreeNodeClass> &BaseKDTree<T, TreeNodeClass>::operator=(
    const BaseKDTree<T, TreeNodeClass> &rhs)
{
    /*std::cerr << "DD: Entering operator=..." << std::endl;;*/
    /* check for self-assignment */
    if (this == &rhs) {
        return *this;
    }
    else {
        copyTree(rhs);
        return *this;
    }
}



template <class T, class TreeNodeClass>
BaseKDTree<T, TreeNodeClass>::BaseKDTree()
{
    setUpEmptyTree();
}





template <class T, class TreeNodeClass>
unsigned int BaseKDTree<T, TreeNodeClass>::size() const
{
    return mySize;
}




template <class T, class TreeNodeClass>
void BaseKDTree<T, TreeNodeClass>::debugPrint() const
{
  std::cout << "KDTree: dims " << myK << std::endl;
  std::vector<double>::iterator myIter;

  std::cout << "UBounds: ";
  for (myIter = myUBounds.begin();
       myIter != myUBounds.end();
       myIter++)
    {
      std::cout << *myIter << " ";
    }

  std::cout << std::endl << "LBounds: ";
  for (myIter = myLBounds.begin();
       myIter != myLBounds.end();
       myIter++)
    {
      std::cout << *myIter << " ";
    }
  std::cout << std::endl;
  myRoot->debugPrint(0);
  
}





template <class T, class TreeNodeClass>
void BaseKDTree<T, TreeNodeClass>::buildFromData(
    std::vector<PointAndValue <T> > 
    pointsAndValues,
    unsigned int k, 
    unsigned int maxLeafSize)
{

    myK = k;

    if (pointsAndValues.size() > 0) 
    {
    
        typename std::vector<PointAndValue<T>,std::allocator<PointAndValue<T> > >::iterator myIter;
        
        std::vector <std::vector<double> > pointsByDimension;
        std::vector <double> pointsUBounds, pointsLBounds;
        double tmpMax, tmpMin;
        std::vector <PointAndValue <T> > pointsAndValuesCopy;
        std::vector <std::vector<double>*> allocatedDoubleVecs;

        /* sanity check */
        if (k < 1) {
            throw LSST_EXCEPT(BadParameterException,
  "EE: KDTree:  k (number of dimensions in data) must be at least 1!\n");
        }
        if (maxLeafSize < 1) {
            throw LSST_EXCEPT(BadParameterException, 
         "EE: KDTree: max leaf size must be strictly positive!\n");
        }
        /* make sure all points from of pointsAndValues are of valid length */
        for (myIter = pointsAndValues.begin(); 
             myIter != pointsAndValues.end(); 
             myIter++) {
            unsigned int pointSize = myIter->getPoint().size();
            if (pointSize < k) {
                throw LSST_EXCEPT(BadParameterException, 
                                  "Got point has size <  k");
            }
        }

        /* find upper and lower bounds of points in each dimension */
  
        for (unsigned int i = 0; i < myK; i++) {
            std::vector<double> *tmpDoubleVec = new std::vector<double>;
            pointsByDimension.push_back(*tmpDoubleVec);
            allocatedDoubleVecs.push_back(tmpDoubleVec);
        }

        for (unsigned int i = 0; i < myK; i++) {

            for (myIter = pointsAndValues.begin();
                 myIter != pointsAndValues.end();
                 myIter++) {
                pointsByDimension[i].push_back(myIter->getPoint()[i]);
            }

            tmpMax = *(std::max_element(pointsByDimension[i].begin(),
                                        pointsByDimension[i].end()));
            tmpMin = *(std::min_element(pointsByDimension[i].begin(),
                                        pointsByDimension[i].end()));
            pointsUBounds.push_back(tmpMax);
            pointsLBounds.push_back(tmpMin);
        }

        for (unsigned int i = 0; i < myK; i++){
            delete allocatedDoubleVecs[i];
        }

        myUBounds = pointsUBounds;
        myLBounds = pointsLBounds;
  
        /* make a copy of the pointsAndValues we were given; it will be
           edited destructively.
        */

        pointsAndValuesCopy = pointsAndValues;
  
        /* create the root of the tree (and the rest of the tree
         * recursively), save it to private var. */
        unsigned int idCounter = 0;
        myRoot = new TreeNodeClass(
            pointsAndValuesCopy, k, maxLeafSize, 0,
            pointsUBounds, pointsLBounds, idCounter);
        // don't set hasData until now, when the tree is actually built.
        mySize = idCounter;
        hasData = true;
    }
}






}} // close namespace lsst::mops

#endif

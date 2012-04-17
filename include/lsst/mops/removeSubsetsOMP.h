// -*- LSST-C++ -*-


/* jonathan myers this is an attempt at doing removal of subset tracklets in a
 * much faster way than the naive algorithm.*/

#ifndef LSST_REMOVE_SUBSETS_H
#define LSST_REMOVE_SUBSETS_H

#include <vector>
#include <set>

#include "lsst/mops/Exceptions.h"

namespace lsst {
namespace mops {


    class SubsetRemover {
    public:

        
        /* jmyers may 17 2011 - 
           
           as per Rob Sinkovits' suggestion, we get our output with a
           vector of equal length to track vector; if outVector[i] ==
           True then trackVector[i] is a track to keep.  This helps
           with speeding up parallelization as well as reducing the
           memory overhead (since such a structure was previously used
           internally to the removeSubsets implementation, to check
           prevent redundant tracks from being written.)

         */
        void removeSubsetsPopulateOutputVector(
            const std::vector<std::set<unsigned int> > *tracksVector, 
            std::vector<bool> &outVector,
            bool shortCircuit=true,
            bool sortBeforeIntersect=false);
        
    };
    
    

}} // close namespace lsst::mops



#endif

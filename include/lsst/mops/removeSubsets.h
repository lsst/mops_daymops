// -*- LSST-C++ -*-


/* jonathan myers this is an attempt at doing removal of subset tracklets in a
 * much faster way than the naive algorithm.*/

#ifndef LSST_REMOVE_SUBSETS_H
#define LSST_REMOVE_SUBSETS_H

#include <vector>

#include "Tracklet.h"


namespace lsst {
namespace mops {


    class SubsetRemover {
    public:
        void removeSubsetsPopulateOutputVector(
            const std::vector<Tracklet> *pairsVector, 
            std::vector<Tracklet> &outVector);
        
    };

    void putLongestPerDetInOutputVector(const std::vector<Tracklet> *pairsVector, 
                                        std::vector<Tracklet> &outputVector);

}} // close namespace lsst::mops



#endif

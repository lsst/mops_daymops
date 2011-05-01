// -*- LSST-C++ -*-


/* jmyers 4/29/2011
 */

/*  
    Linkage is the base class for Track and Tracklet; it should export
    the same interface, allowing the same Tree class to hold Tracks or
    Tracklets.

    You should be able to differentiate between them easily by
    checking the dimensionality of their best-fit functions; tracklets
    have linear underlying models and tracks have quadratic or higher.
 */

#ifndef LSST_LINKAGEVECTOR_H
#define LSST_LINKAGEVECTOR_H

#include "lsst/mops/Linkage.h"

namespace lsst {
namespace mops {

class LinkageVector {
public:
     virtual Linkage* at(unsigned int i) 
     {
	  panic();
	  return NULL;
     }

    virtual unsigned int size() const
        {
            panic();
            return 0;
        }

    virtual void push_back(Linkage* l)
        {
            panic();
        }


private:
     void panic() const { throw LSST_EXCEPT(InheritanceException,
      "virtual base class function was called in LinkageVector!\n");
     }

};


} } //close namespace lsst::mops

#endif

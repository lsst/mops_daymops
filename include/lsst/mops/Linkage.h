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

#ifndef LSST_LINKAGE_H
#define LSST_LINKAGE_H

#include <vector>
#include <set>
#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/Exceptions.h"


namespace lsst {
namespace mops {

class Linkage {
public:
    
    /* Linkages should demand calculateFitFunction be called before
     * getFitFunction. */
    virtual void calculateFitFunction() 
        { panic(); return; }
    
    /* functions should be modified to hold constants of the motions
       of equation going from low-order to high, e.g.  [p0, vel] or [p0,
       vel] */
    virtual void getFitFunction(double &epoch,
                                std::vector<double> &raFunc,
                                std::vector<double> &fitFunc) 
        { panic(); return; }
    
    virtual const std::set<unsigned int> getComponentDetectionIndices() const 
        { panic(); std::set<unsigned int> dummy; return dummy;}
    
    
    virtual double getStartTime(
        const std::vector<MopsDetection> &dets) const
        { panic(); return 0.; }
    
    virtual unsigned int getId() const { panic(); return 0; };
    virtual void setId(unsigned int newId) { panic(); };

    virtual double getDeltaTime(
        const std::vector<MopsDetection> &dets) const
        { panic(); return 0.;}

private:
     void panic() const { throw LSST_EXCEPT(InheritanceException,
      "virtual base class function was called in Linkage!\n");
     }

};


} } //close namespace lsst::mops

#endif

// -*- LSST-C++ -*-


/* jmyers 2/11/11
 */


#include "lsst/mops/Tracklet.h"

namespace lsst {
namespace mops {

Tracklet::Tracklet()
{
     isCollapsed = false; 
     velocityRA = 0.; 
     velocityDec = 0.;
}



Tracklet::Tracklet(std::set <unsigned int> startIndices) 
{ 
     isCollapsed = false; 
     indices = startIndices;
}
     


double Tracklet::getStartTime(std::vector<MopsDetection> dets) 
{
     if (indices.size() < 1) {
	  LSST_EXCEPT(UninitializedException, 
		      "Tracklet: start time requested, but tracklet contains no detections.");
     }
          
     bool haveFirst = false;
     double firstTime = 0.;

     std::set<unsigned int>::const_iterator indIt;
     for (indIt = indices.begin(); indIt != indices.end(); indIt++) {

	  if (*indIt > dets.size()) {
	       LSST_EXCEPT(BadIndexException, 
  "Tracklet: the detection vector we were handed cannot possibly go with this tracklet.");
	  }

	  double thisTime = dets.at(*indIt).getEpochMJD();
	  if ((!haveFirst) || (thisTime < firstTime)) {
	       firstTime = thisTime;
	       haveFirst = true;
	  }
     }
     return firstTime;
}

}}

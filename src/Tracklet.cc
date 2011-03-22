// -*- LSST-C++ -*-


/* jmyers 2/11/11
 */


#include "lsst/mops/Tracklet.h"

namespace lsst {
namespace mops {

Tracklet::Tracklet()
{
    myId = -1;
    isCollapsed = false; 
    velocityRA = 0.; 
    velocityDec = 0.;
}



Tracklet::Tracklet(std::set <unsigned int> startIndices) 
{ 
    myId = -1;
    velocityRA = 0.;
    velocityDec = 0.;
    isCollapsed = false; 
    indices = startIndices;
}
     

MopsDetection Tracklet::getFirstDetection(
    const std::vector<MopsDetection> &dets) const
{
     if (indices.size() < 1) {
	  LSST_EXCEPT(UninitializedException, 
   "Tracklet: start time requested, but tracklet contains no detections.");
     }
          
     bool haveFirst = false;
     MopsDetection firstDet;
     double firstTime;

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
               firstDet = dets.at(*indIt);
	  }
     }
     return firstDet;    
}



double Tracklet::getStartTime(const std::vector<MopsDetection> &dets) const
{
    return getFirstDetection(dets).getEpochMJD();
}



double Tracklet::getDeltaTime(const std::vector<MopsDetection> &dets) const
{
    return getLastDetection(dets).getEpochMJD() - getStartTime(dets);
}


MopsDetection Tracklet::getLastDetection(
    const std::vector<MopsDetection> &dets)
    const
{
     if (indices.size() < 1) {
	  LSST_EXCEPT(UninitializedException, 
   "Tracklet: start time requested, but tracklet contains no detections.");
     }
          
     bool haveLast = false;
     MopsDetection lastDet;
     double lastTime;

     std::set<unsigned int>::const_iterator indIt;
     for (indIt = indices.begin(); indIt != indices.end(); indIt++) {

	  if (*indIt > dets.size()) {
	       LSST_EXCEPT(BadIndexException, 
  "Tracklet: the detection vector we were handed cannot possibly go with this tracklet.");
	  }

	  double thisTime = dets.at(*indIt).getEpochMJD();
	  if ((!haveLast) || (thisTime > lastTime)) {
	       lastTime = thisTime;
	       haveLast = true;
               lastDet = dets.at(*indIt);
	  }
     }
     return lastDet;    

}


Tracklet & Tracklet::operator= (const Tracklet &other)
 {
    copy(other);
    return *this;
 } 

bool Tracklet::operator==(const Tracklet &other) const
{
    if (indices == other.indices){
        return true;
    }
    else {
        return false;
    }
}

bool Tracklet::operator!=(const Tracklet &other) const
{
    if (indices != other.indices){
        return true;
    }
    else {
        return false;
    }
    
}

bool Tracklet::operator<(const Tracklet &other) const
{
    return (indices < other.indices);
}

void Tracklet::copy(const Tracklet &other)
{
    indices = other.indices;
    velocityRA = other.velocityRA;
    velocityDec = other.velocityDec;
    isCollapsed = other.isCollapsed;
    myId = other.myId;
}

}}

// -*- LSST-C++ -*-
/* jonathan myers 
   4/08/10
*/

#include <set>


#include "lsst/mops/TrackletVector.h"



namespace lsst { 
namespace mops {




void TrackletVector::push_back(const Tracklet &newTracklet) {
     componentTracklets.push_back(newTracklet);
}


Tracklet TrackletVector::at(unsigned int i) const {
    return componentTracklets.at(i);
}



unsigned int TrackletVector::size() const {
     return componentTracklets.size();
}


    
bool TrackletVector::isSubsetOf(const TrackletVector &other) const {
     if (other.size() < this->size()) {
	  return false;
     }
     else {
         // create a copy TrackletVector as a set for fast searching
         std::set<Tracklet> otherTracklets;
         std::vector<Tracklet>::const_iterator copyIter;
         for (copyIter = other.componentTracklets.begin();
              copyIter != other.componentTracklets.end();
              copyIter++) {
             otherTracklets.insert(*copyIter);
         }
         
         // now see if the set contains all of our tracklets
         std::vector<Tracklet>::const_iterator tIter;
         for (tIter = componentTracklets.begin();
              tIter != componentTracklets.end();
              tIter++) {
             if (otherTracklets.find(*tIter) 
                 == otherTracklets.end()) {
                 return false;
             }
         }
         return true;
     }
}



bool TrackletVector::operator==(const TrackletVector other) const {

    // if two sets are subsets of each other, they are equal.

    return (isSubsetOf(other) && other.isSubsetOf(*this));

}

bool TrackletVector::operator!=(const TrackletVector &other) const {
     return ! (*this == other);
}



}} // close namespace lsst::mops

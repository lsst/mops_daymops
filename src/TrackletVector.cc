// -*- LSST-C++ -*-
/* jonathan myers 
   4/08/10
*/

#include <ios>
#include <set>

#include "lsst/mops/Exceptions.h"
#include "lsst/mops/fileUtils.h"
#include "lsst/mops/TrackletVector.h"



namespace lsst { 
namespace mops {



TrackletVector::TrackletVector(std::string outFileName, bool useCache, unsigned int cacheSize) 
{
    if (useCache) {
        this->useCache = true;
        this->cacheSize = cacheSize;
    }
    else {
        this->useCache = false;
        this->cacheSize = 0;        
    }

    useOutFile = true;
    outFile.open(outFileName.c_str(), std::ios_base::out | std::ios_base::app);

}





void TrackletVector::purgeToFile() 
{
    if (useOutFile) {
        if (componentTracklets.size() != 0) {
            writeToFile();
        }        
    }
    else {
        throw LSST_EXCEPT(BadParameterException,
                          "TrackletVector: Cannot purge to file in a trackletVector created without a file.");
    }

    if (useCache) {
        componentTracklets.clear();
    }
}






TrackletVector::~TrackletVector() 
{
    if (useOutFile) {
        purgeToFile();
        outFile.close();
    }
}





void TrackletVector::writeToFile()
{
    // TBD: eventually, get rid of fileUtils entirely and put its code here.
    writeTrackletsToOutFile(&componentTracklets, outFile);
}





void TrackletVector::push_back(const Tracklet &newTracklet) {


    if ((useOutFile) && (useCache) && (componentTracklets.size() >= cacheSize)) {
        purgeToFile();
    }

    componentTracklets.push_back(newTracklet);
}






Tracklet TrackletVector::at(unsigned int i) const {

    if (useCache) {
      throw LSST_EXCEPT(KnownShortcomingException,  
                          "trackletVector: Cannot call do random access via at() when using output cacheing.");
    }

    return componentTracklets.at(i);
}





unsigned int TrackletVector::size() const {

    if (useCache) {
        throw LSST_EXCEPT(KnownShortcomingException, 
                          "trackletVector: Cannot request 'size' when using cacheing.");
    }

    return componentTracklets.size();
}



    
bool TrackletVector::isSubsetOf(const TrackletVector &other) const {

    if (useCache) {
      throw LSST_EXCEPT(KnownShortcomingException, 
                          "trackletVector: Cannot call isSubsetOf when using cacheing.");
    }

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



bool TrackletVector::operator==(const TrackletVector &other) const {

    if (useCache) {
        throw LSST_EXCEPT(KnownShortcomingException, 
                          "trackletVector: Cannot call == when using cacheing.");
    }

    // if two sets are subsets of each other, they are equal.

    return (isSubsetOf(other) && other.isSubsetOf(*this));

}

bool TrackletVector::operator!=(const TrackletVector &other) const {

    if (useCache) {
        throw LSST_EXCEPT(KnownShortcomingException, 
                          "trackletVector: Cannot call != when using cacheing.");
    }

    return ! (*this == other);
}



}} // close namespace lsst::mops

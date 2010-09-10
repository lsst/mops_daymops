// -*- LSST-C++ -*-
/* jonathan myers 
   4/08/10
*/


#include "lsst/mops/Exceptions.h"
#include "lsst/mops/TrackSet.h"

namespace lsst { 
namespace mops {




TrackSet::TrackSet(std::string outFileName, bool useCache, unsigned int cacheSize) 
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

void TrackSet::purgeToFile() 
{
    if (useOutFile) {
        if (componentTracks.size() != 0) {
            writeToFile();
        }        
    }
    else {
        throw LSST_EXCEPT(BadParameterException,
                          "TrackSet: Cannot purge to file in a trackSet created without a file.");
    }

    if (useCache) {
        componentTracks.clear();
    }
}






TrackSet::~TrackSet() 
{
    if (useOutFile) {
        purgeToFile();
        outFile.close();
    }
}





void TrackSet::writeToFile()
{
    std::set<Track>::const_iterator curTrack;
    for (curTrack = componentTracks.begin(); 
         curTrack != componentTracks.end();
         curTrack++) {
        std::set<unsigned int>::const_iterator detIter;
	
        for (detIter = curTrack->componentDetectionIndices.begin();
             detIter != curTrack->componentDetectionIndices.end();
             detIter++) {
            outFile << *detIter << " ";
        }
        outFile << std::endl;
    }
    outFile.close();
}






void TrackSet::insert(const Track &newTrack) {
    componentTracks.insert(newTrack);
}



unsigned int TrackSet::size() const {
    return componentTracks.size();
}



bool TrackSet::isSubsetOf(const TrackSet &other) const {
    if (other.size() < this->size()) {
        return false;
    }
    else {
        // create a copy TrackSet
        std::set<Track>::const_iterator tIter;
        for (tIter = componentTracks.begin();
             tIter != componentTracks.end();
             tIter++) {
            if (other.componentTracks.find(*tIter) 
                == other.componentTracks.end()) {
                return false;
            }
        }
        return true;
    }
}



bool TrackSet::operator==(const TrackSet &other) const {
    return (componentTracks == other.componentTracks);
}



bool TrackSet::operator!=(const TrackSet &other) const {
    return ! (*this == other);
}

}} // close namespace lsst::mops

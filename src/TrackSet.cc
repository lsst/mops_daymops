// -*- LSST-C++ -*-
/* jonathan myers 
   4/08/10
*/


#include "lsst/mops/Exceptions.h"
#include "lsst/mops/TrackSet.h"

#include <iomanip>

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
        std::set<unsigned int> diaIds = curTrack->getComponentDetectionDiaIds();
        for (detIter = diaIds.begin();
             detIter != diaIds.end();
             detIter++) {
            outFile << *detIter << " ";
        }
        /* jmyers: new for ticket 1415: also append best-fit function data.
           note that there is currently no mechanism for checking that this
           information has been calculated is or up-to-date; this is OK for now
           because this routine is only ever called by linkTracklets, which
           calculates an up-to-date fit function for every track, after adding
           support points.
         */
        double epoch, ra0, raV, raAcc, dec0, decV, decAcc;
        (*curTrack).getBestFitQuadratic(epoch, ra0, raV, raAcc, 
                                     dec0, decV, decAcc);
        // outFile << std::setprecision(12) 
        //         << std::scientific 
        //         << "epoch=" << epoch 
        //         << " ra0=" << ra0 << " raV=" << raV << " raAcc=" << raAcc
        //         << " dec0=" << dec0 << " decV=" << decV << " decAcc=" << decAcc;
        
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

// -*- LSST-C++ -*-
/* jonathan myers 
   4/08/10
*/

#ifndef LSST_TRACKSET_H
#define LSST_TRACKSET_H

#include <iostream>
#include <fstream>
#include <set>

#include "Track.h"

namespace lsst {
namespace mops {


class TrackSet {
public:
    // create a normal TrackSet, which is just a container class. no file behaviors.
    TrackSet() { useCache = false; cacheSize = 0; useOutFile = false;};

    /*
     * if useCache == True: create a TrackSet which holds at most cacheSize 
     * elements, and if elements > that number are added, then write 
     * them to file with name outFileName.
     *
     * if useCache == False, open file with name outFileName for writing.
     * on destruction or call to purgeCacheToFile, write all contents to that outfile.
     */
    TrackSet(std::string outFileName, bool useCache=false, unsigned int cacheSize=0);


    /*
     * calls purgeToFile for you, in case there are still items in the
     * cache that need to be written to disk.  Closes the outFile if in use.
     */

    ~TrackSet();

    /*
     * if there is no outfile, raise exception.
     *
     * if there is an outfile, write all contents to that outfile.
     * 
     * if cacheing is enabled, clear the local entries as well.
     */
    void purgeToFile();

    std::set<Track> componentTracks;

    void insert(const Track &newTrack);

    unsigned int size() const;
    
    /* These operators will fail with exception if useCache is true. They are
     * useful for debugging and unit testing. */

    bool isSubsetOf(const TrackSet &other) const;
    
    bool operator==(const TrackSet &other) const;
    
    bool operator!=(const TrackSet &other) const;

private:
    void writeToFile();
    bool useCache;
    std::ofstream outFile;
    bool useOutFile;
    unsigned int cacheSize;
    

};



}} // close namespace lsst::mops

#endif

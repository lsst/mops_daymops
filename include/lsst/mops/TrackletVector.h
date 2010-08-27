// -*- LSST-C++ -*-
/* jonathan myers 
   4/08/10
*/

#ifndef LSST_TRACKLETVECTOR_H
#define LSST_TRACKLETVECTOR_H

#include <fstream>
#include <string>
#include <vector>


#include "Tracklet.h"

namespace lsst {
namespace mops {



// formerly inherited from Persistable, but for now we aren't using
// that functionality.

class TrackletVector {
public:

    // create a normal TrackletVector, which is just a container class. no file behaviors.
    TrackletVector() { useCache = false; cacheSize = 0; useOutFile = false;};

    /*
     * if useCache == True: create a TrackletVector which holds at most cacheSize 
     * elements, and if elements > that number are added, then write 
     * them to file with name outFileName.
     *
     * if useCache == False, open file with name outFileName for writing.
     * on destruction or call to purgeCacheToFile, write all contents to that outfile.
     */
    TrackletVector(std::string outFileName, bool useCache=false, unsigned int cacheSize=0);


    /*
     * calls purgeToFile for you, in case there are still items in the
     * cache that need to be written to disk.  Closes the outFile if in use.
     */

    ~TrackletVector();

    /*
     * if there is no outfile, raise exception.
     *
     * if there is an outfile, write all contents to that outfile.
     * 
     * if cacheing is enabled, clear the local entries as well.
     */
    void purgeToFile();

    void push_back(const Tracklet &newTracklet);

    Tracklet at(unsigned int) const;

    unsigned int size() const;
    


    /*
     * see warnings in Tracklet.h.  If two tracklet sets were generated from
     * different detection vectors, even if those vectors contain the same
     * detections in a different ordering, the results of these comparison
     * functions are NOT meaningful!
     */

    /*
      this will return true iff this the other vector contains every tracklet
      that we do.  ordering is ignored.  

      to make the implementation FAST, we will create a COPY of the other
      trackletVector's tracklets represented as a set; be wary of memory usage.
     */
    bool isSubsetOf(const TrackletVector &other) const;


    /*
     * two vectors are == if and only if they contain the same tracklets;
     * ordering is NOT considered.  For fast performance, all Tracklets in other
     * will be copied, then deleted, then all tracklets in this TrackletVector
     * will be copied, then deleted.
     */
    bool operator==(const TrackletVector &other) const;

    bool operator!=(const TrackletVector &other) const;

private:
    void writeToFile();
    std::vector<Tracklet> componentTracklets;
    bool useCache;
    std::ofstream outFile;
    bool useOutFile;
    unsigned int cacheSize;
    

};

}} // close namespace lsst::mops



#endif

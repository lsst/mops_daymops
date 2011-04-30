// -*- LSST-C++ -*-

#include "lsst/mops/daymops/linkTracklets/linkTracklets.h"
#include "lsst/mops/Track.h"
#include "lsst/mops/TrackVector.h"
#include "lsst/mops/fileUtils.h"

#include <string>
#include <iostream>


namespace lsst { namespace mops {



int readTracksWriteStats(int argc, char** argv) {

    
    if (argc != 4) {
        std::string usage("USAGE: readTracksWriteStats <inDets> <inTracks_byIndices> <outFile> \n");
        std::cout << usage;
        std::cerr << usage;
        return 1;
    }
    
    std::string inDetsFile(argv[1]);
    std::string inFile(argv[2]);
    std::string outFile(argv[3]);
    


    // need to set observatory location
    linkTrackletsConfig searchConfig;
    lsst::mops::MopsDetection::setObservatoryLocation(searchConfig.obsLat, searchConfig.obsLong);
     const double astromErr =  searchConfig.defaultAstromErr;
     std::cerr << "Using defaultAstromErr " << astromErr << '\n';


     
     TrackVector myTv;
     std::vector<MopsDetection> allDets;
     std::cout << "Reading dets...\n";
     populateDetVectorFromFile(inDetsFile, allDets, astromErr);
     std::cout << "Reading input tracks (as indices)...\n";
     myTv.populateFromFile(inFile, allDets);
     std::cout << "Calculating topocentric correction of detections...\n";
     calculateTopoCorr(allDets, searchConfig);
     std::cout << "Calculating stats, writing output tracks (as Dias) and stats...\n";
     myTv.writeTracksAndStatsToFile(outFile, allDets);
     std::cout << "Done.\n";

     return 0;
}

}} // close lsst::mops




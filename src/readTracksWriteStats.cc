// -*- LSST-C++ -*-

#include "lsst/mops/Track.h"
#include "lsst/mops/TrackVector.h"
#include "lsst/mops/fileUtils.h"
#include "lsst/mops/daymops/linkTracklets/linkTracklets.h"

#include <string>
#include <iostream>


namespace lsst { namespace mops {

int doIt(int argc, char** argv) {


     if (argc != 4) {
	  std::string usage("USAGE: readTracksWriteStats <inDets> <inTracks_byIndices> <outFile> \n");
	  std::cout << usage;
	  std::cerr << usage;
	  return 1;
     }

     lsst::mops::linkTrackletsConfig searchConfig;

     std::string inDetsFile(argv[1]);
     std::string inFile(argv[2]);
     std::string outFile(argv[3]);
     
     TrackVector myTv;
     std::vector<MopsDetection> allDets;
     std::cout << "Reading dets...\n";
     populateDetVectorFromFile(inDetsFile, allDets, searchConfig.defaultAstromErr);
     for (uint i = 0; i < allDets.size(); i++) {
         allDets[i].setObservatoryLocation(searchConfig.obsLat, searchConfig.obsLong);
         allDets[i].calculateTopoCorr();
     }
     std::cout << "First det has id, ra, dec " << allDets[0].getID() << " "
               <<  allDets[0].getRA() << " " << allDets[0].getDec() << "\n";
     std::cout << "Reading input tracks (as indices)...\n";
     myTv.populateFromFile(inFile, allDets);
     std::cout << "Writing output tracks (as Dias) and stats...\n";

     myTv.writeTracksToFile(outFile);

     myTv.writeTracksAndStatsToFile(outFile, allDets);
     std::cout << "Done.\n";

     return 0;
}

}} // close lsst::mops




int main(int argc, char** argv) 
{
     return lsst::mops::doIt(argc, argv);

}

// -*- LSST-C++ -*-

#include "lsst/mops/Track.h"
#include "lsst/mops/TrackVector.h"
#include "lsst/mops/fileUtils.h"

#include <string>
#include <iostream>


namespace lsst { namespace mops {

int doIt(int argc, char** argv) {

     /* compilation fails unless there is an instance of Track in this
      * program. I have no idea why.  Without it, g++ will claim that
      * none of Track's member functions are found, and point to their
      * calls in TrackVector. ugh!  */
     Track unusedMagicTrack; 


     if (argc != 4) {
	  std::string usage("USAGE: readTracksWriteStats <inDets> <inTracks_byIndices> <outFile> \n");
	  std::cout << usage;
	  std::cerr << usage;
	  return 1;
     }

     std::string inDetsFile(argv[1]);
     std::string inFile(argv[2]);
     std::string outFile(argv[3]);
     
     TrackVector myTv;
     std::vector<MopsDetection> allDets;
     std::cout << "Reading dets...\n";
     populateDetVectorFromFile(inDetsFile, allDets);
     std::cout << "Reading input tracks (as indices)...\n";
     myTv.populateFromFile(inFile, allDets);
     std::cout << "Writing output tracks (as Dias) and stats...\n";
     myTv.writeTracksAndStatsToFile(outFile, allDets);
     std::cout << "Done.\n";

     return 0;
}

}} // close lsst::mops




int main(int argc, char** argv) 
{
     return lsst::mops::doIt(argc, argv);

}

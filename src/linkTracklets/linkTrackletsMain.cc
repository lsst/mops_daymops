#include <boost/lexical_cast.hpp>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <ctime>
#include <time.h>

#include "lsst/mops/daymops/linkTracklets/linkTracklets.h"
#include "lsst/mops/fileUtils.h"




namespace lsst {
     namespace mops {




}} // close lsst::mops

int main(int argc, char* argv[])
{

     double dif;
     time_t start = time(NULL);

     lsst::mops::linkTrackletsConfig searchConfig; 

     /* 
      * we know we're being run from the command line, so set verbosity high.
      */

     searchConfig.myVerbosity.printStatus = true;
     searchConfig.myVerbosity.printVisitCounts = true;
     searchConfig.myVerbosity.printTimesByCategory = true;
     searchConfig.myVerbosity.printBoundsInfo = true;
     
     int bufferSize = 1000;

     std::string helpString = 
	  std::string("Usage: linkTracklets -d <detections file> -t <tracklets file> -o <output (tracks) file>") + std::string("\n") +
	  std::string("  optional arguments: ") + std::string("\n") +
	  std::string("     -e / --detectionErrorThresh (float) : maximum allowed observational error, default = ")
	  + boost::lexical_cast<std::string>(searchConfig.detectionLocationErrorThresh) + std::string("\n") +
	  std::string("     -D / --maxDecAcceleration (float) : maximum sky-plane acceleration of a track (declination),  default = ")
	  + boost::lexical_cast<std::string>(searchConfig.maxDecAccel) + std::string("\n") +
	  std::string("     -R / --maxRAAcceleration (float) : maximum sky-plane acceleration of a track (RA), default = ")
	  + boost::lexical_cast<std::string>(searchConfig.maxRAAccel) +  std::string("\n") +
	  std::string("     -F / --latestFirstEndpoint (float) : if specified, only search for tracks with first endpoint before time specified")
	  + std::string("\n") +
	  std::string("     -L / --earliestLastEndpoint (float) : if specified, only search for tracks with last endpoint after time specified")
	  +  std::string("\n") +
	  std::string("     -u / --minNights (int) : require tracks contain detections from at least this many nights, default = ")
	  + boost::lexical_cast<std::string>(searchConfig.minUniqueNights) +  std::string("\n") +
	  std::string("     -s / --minDetections (int) : require tracks contain at least this many detections, default = ")
	  + boost::lexical_cast<std::string>(searchConfig.minDetectionsPerTrack) +  std::string("\n") +
	  std::string("     -b / --outputBufferSize (int) : number of tracks to buffer in memory before flushing output. Default = ")
	  + boost::lexical_cast<std::string>(bufferSize) +  std::string("\n") +
	  std::string("     -n / --leafNodeSize (int) : set max leaf node size for nodes in KDTree")
	  +  std::string("\n");

     static const struct option longOpts[] = {
	  { "detectionsFile", required_argument, NULL, 'd' },
	  { "trackletsFile", required_argument, NULL, 't' },
	  { "outputFile", required_argument, NULL, 'o' },
	  { "detectionErrorThresh", required_argument, NULL, 'e'},
	  { "maxDecAcceleration", required_argument, NULL, 'D'},
	  { "maxRAAcceleration", required_argument, NULL, 'R'},
	  { "latestFirstEndpoint", required_argument, NULL, 'F'},
	  { "earliestLastEndpointTime", required_argument, NULL, 'L'},
	  { "minNights", required_argument, NULL, 'u'},
	  { "minDetections", required_argument, NULL, 's'},
	  { "outputBufferSize", required_argument, NULL, 'b'},
	  { "leafNodeSize", required_argument, NULL, 'n'},
	  { "help", no_argument, NULL, 'h' },
	  { NULL, no_argument, NULL, 0 }
     };  
     
     
     std::stringstream ss;
     std::string detectionsFileName = "";
     std::string trackletsFileName = "";

     
     int longIndex = -1;
     const char *optString = "d:t:o:e:D:R:F:L:u:s:b:n:h";
     int opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
     while( opt != -1 ) {
	  switch( opt ) {
	  case 'd':	       
	       detectionsFileName = optarg;
	       break;
	  case 't':
	       trackletsFileName = optarg;
	       break;
	  case 'o':
	       searchConfig.outputFile = optarg;
	       break;
	  case 'e':
	       searchConfig.detectionLocationErrorThresh = atof(optarg);
	       break;
	  case 'D':
	       searchConfig.maxDecAccel = atof(optarg);

	       break;
	  case 'R':
	       searchConfig.maxRAAccel = atof(optarg);
	       break;
	  case 'F':
	       searchConfig.restrictTrackStartTimes = true;
	       searchConfig.latestFirstEndpointTime = atof(optarg);
	       std::cout << "Got latest first endpoint time = " << 
		    std::setprecision(12) << searchConfig.latestFirstEndpointTime
			 << std::endl;
	       break;
	  case 'L':
	       searchConfig.restrictTrackEndTimes = true;
	       searchConfig.earliestLastEndpointTime = atof(optarg);
	       break;
	  case 'u':
	       searchConfig.minUniqueNights = atoi(optarg);
	       std::cout << "Set min unique nights per track: " 
			 << searchConfig.minUniqueNights << "\n";
	       break;
	  case 's':
	       searchConfig.minDetectionsPerTrack = atoi(optarg);
	       std::cout << "Set min detections per track: " 
			 << searchConfig.minDetectionsPerTrack << "\n";
	       break;
	  case 'b':
	       bufferSize = atoi(optarg);
	       std::cout << " Set output buffer size = " 
			 << bufferSize << std::endl;
	       if (bufferSize < 1) {
		    std::cerr << "Illegal output buffer size. Exiting.\n";
		    return -1;
	       }
	       break;
	  case 'n':
	       searchConfig.leafSize = atoi(optarg);
	       std::cout << " Set leaf node size = " 
			 << searchConfig.leafSize << std::endl;
	       break;
	  case 'h':
	       std::cout << helpString << std::endl;
	       return 0;
	  default:
	       break;
	  }
	  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
     }

// Set static obsLat and obsLong in MopsDetection

     lsst::mops::MopsDetection::setObservatoryLocation(searchConfig.obsLat, searchConfig.obsLong);

     if ((detectionsFileName == "") || (trackletsFileName == "")) {
	  std::cerr << helpString << std::endl;
	  return 1;
     }

     std::vector<lsst::mops::MopsDetection> allDets;
     std::vector<lsst::mops::Tracklet> allTracklets;
     lsst::mops::TrackSet * resultTracks;
     searchConfig.outputMethod = lsst::mops::IDS_FILE_WITH_CACHE;
     searchConfig.outputBufferSize = bufferSize;
     
     
     const double astromErr =  searchConfig.defaultAstromErr;
     std::cerr << "Using defaultAstromErr " << astromErr << '\n';
     std::cout << "Reading detections file " << std::endl;
     populateDetVectorFromFile(detectionsFileName, allDets, astromErr);
     std::cout << "Calculating topocentric correction " << std::endl;
     calculateTopoCorr(allDets, searchConfig);
     std::cout << "Reading tracklets file " << std::endl;
     populatePairsVectorFromFile(trackletsFileName, allTracklets);

     dif = lsst::mops::timeElapsed(start);
     std::cout << "Reading input took " << std::fixed << std::setprecision(10) 
	       <<  dif  << " seconds." <<std::endl;     

     

     resultTracks = lsst::mops::linkTracklets(allDets, allTracklets, searchConfig);

     
     resultTracks->purgeToFile();
     std::cout << "Results successfully written to disk." << std::endl;
     
     std::cout << "Done. Exiting successfully." << std::endl;

     dif = lsst::mops::timeElapsed(start);
     std::cout << "Completed after " << std::fixed << std::setprecision(10) 
	       <<  dif  << " seconds." <<std::endl;     
     lsst::mops::printMemUse();
     return 0;

    	    
}

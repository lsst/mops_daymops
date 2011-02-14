#include <boost/lexical_cast.hpp>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <getopt.h>
#include <iomanip>


#include "lsst/mops/daymops/linkTracklets/linkTracklets.h"
#include "lsst/mops/fileUtils.h"



#define PRINT_TIMING_INFO true

#ifdef PRINT_TIMING_INFO
#include <ctime>
#endif

namespace lsst {
     namespace mops {

double timeElapsed(clock_t priorEvent)
{
     return ( std::clock() - priorEvent ) / (double)CLOCKS_PER_SEC;
}



}} // close lsst::mops

int main(int argc, char* argv[])
{

     lsst::mops::linkTrackletsConfig searchConfig; 

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
	  { "help", no_argument, NULL, 'h' },
	  { NULL, no_argument, NULL, 0 }
     };  
     
     
     std::stringstream ss;
     std::string detectionsFileName = "";
     std::string trackletsFileName = "";

     
     int longIndex = -1;
     const char *optString = "d:t:o:e:D:R:F:L:h";
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
	       break;
	  case 'L':
	       searchConfig.restrictTrackEndTimes = true;
	       searchConfig.earliestLastEndpointTime = atof(optarg);
	       break;
	  case 'h':
	       std::cout << helpString << std::endl;
	       return 0;
	  default:
	       break;
	  }
	  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
     }

     if ((detectionsFileName == "") || (trackletsFileName == "")) {
	  std::cerr << helpString << std::endl;
	  return 1;
     }

     std::vector<lsst::mops::MopsDetection> allDets;
     std::vector<lsst::mops::Tracklet> allTracklets;
     lsst::mops::TrackSet * resultTracks;
     searchConfig.outputMethod = lsst::mops::IDS_FILE_WITH_CACHE;
     searchConfig.outputBufferSize = 1000;

     clock_t last;
     double dif;
     if(PRINT_TIMING_INFO) {     
	  last = std::clock();
     }
     populateDetVectorFromFile(detectionsFileName, allDets);
     populatePairsVectorFromFile(trackletsFileName, allTracklets);

     if(PRINT_TIMING_INFO) {     	  
	  dif = lsst::mops::timeElapsed(last);
	  std::cout << "Reading input took " << std::fixed << std::setprecision(10) 
		    <<  dif  << " seconds." <<std::endl;     
     }

     
     if(PRINT_TIMING_INFO) {     
	  last = std::clock();
     }

     resultTracks = lsst::mops::linkTracklets(allDets, allTracklets, searchConfig);

     if(PRINT_TIMING_INFO) {     
	  dif = lsst::mops::timeElapsed (last);
	  std::cout << "linking took " << std::fixed << std::setprecision(10) <<  dif 
		    << " seconds."<<std::endl;     
     }

     if(PRINT_TIMING_INFO) {     
	  last = std::clock();
     }

     resultTracks->purgeToFile();
     std::cout << "Results successfully written to disk." << std::endl;
     

     if(PRINT_TIMING_INFO) {     	  
	  dif = lsst::mops::timeElapsed(last);
	  std::cout << "Writing output took " << std::fixed << std::setprecision(10) 
		    <<  dif  << " seconds." <<std::endl;     
     }

     std::cout << "Done. Exiting successfully." << std::endl;


     return 0;

    	    
}

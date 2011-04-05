// -*- LSST-C++ -*-
#include "mpi.h"
#include <boost/lexical_cast.hpp>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <getopt.h>
#include <iomanip>


//#include "linkTracklets.h"
#include "lsst/mops/rmsLineFit.h"
#include "lsst/mops/daymops/linkTracklets/linkTracklets.h"
#include "lsst/mops/Exceptions.h"
#include "lsst/mops/KDTree.h"
#include "lsst/mops/daymops/linkTracklets/TrackletTree.h"
#include "../../include/lsst/mops/fileUtils.h"


#define PRINT_TIMING_INFO false

#ifdef PRINT_TIMING_INFO
#include <ctime>
#endif

MPI_Datatype bruteForceArgs;


double timeElapsed(clock_t priorEvent)
{
     return ( std::clock() - priorEvent ) / (double)CLOCKS_PER_SEC;
}



int main(int argc, char* argv[])
{
  
  /* 
   * Initialize MPI runtime environment
   */
  int rc = MPI_Init(&argc, &argv);
  if( rc != MPI_SUCCESS ){
    std::cerr << "MPI failed to initialize, aborting." << std::endl;
    exit(rc);
  }

  /*
   * Establish MPI-related variable values
   */
  int rank, numProcessors;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcessors); //number of processors
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); //my rank

  //master node prepares all sets for processing and distributes the workload
  //to the slave nodes, which are all waiting in doLinkingRecurse2
  //processor one is the only one to read data and process arguments
  
  std::vector<lsst::mops::MopsDetection> allDets;
  std::vector<lsst::mops::Tracklet> allTracklets;
  //std::vector<Track> resultTracks;
  lsst::mops::TrackSet resultTracks;
  clock_t last;
  double dif;
  std::string outputFileName = "";
  lsst::mops::linkTrackletsConfig searchConfig; 

  searchConfig.myVerbosity.printStatus = true;
  searchConfig.myVerbosity.printVisitCounts = true;
  searchConfig.myVerbosity.printTimesByCategory = true;
  searchConfig.myVerbosity.printBoundsInfo = true;
  
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
	  std::string("     -n / --leafNodeSize (int) : set max leaf node size for nodes in KDTree")
	  +  std::string("\n");

  
  static const struct option longOpts[] = {
    { "detectionsFile", required_argument, NULL, 'd' },
    { "trackletsFile", required_argument, NULL, 't' },
    { "outputFile", required_argument, NULL, 'o' },
    { "detectionErrorThresh", required_argument, NULL, 'e'},
    { "velocityErrorThresh", required_argument, NULL, 'v'},
    { "maxDecAcceleration", required_argument, NULL, 'D'},
    { "maxRAAcceleration", required_argument, NULL, 'R'},
    { "latestFirstEndpoint", required_argument, NULL, 'F'},
    { "earliestLastEndpointTime", required_argument, NULL, 'L'},
    { "leafNodeSize", required_argument, NULL, 'n'},
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
  };  
  
  
  std::stringstream ss;
  std::string detectionsFileName = "";
  std::string trackletsFileName = "";
  
  
  int longIndex = -1;
  //const char *optString = "d:t:o:e:v:D:R:h";
  const char *optString = "d:t:o:e:D:R:F:L:h:v:n:";
  int opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  while( opt != -1 ) {
    switch( opt ) {
    case 'd':
      detectionsFileName = optarg;
      break;
    case 't':
      trackletsFileName = optarg;
      break;
      /*case 'o':
	ss << optarg;
	ss >> outputFileName;/
	outputFileName = optarg;
	break;*/
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
        std::cerr << "Got latest first endpoint time = " << 
            std::setprecision(12) << searchConfig.latestFirstEndpointTime
                  << std::endl;
        break;
    case 'L':
        searchConfig.restrictTrackEndTimes = true;
        searchConfig.earliestLastEndpointTime = atof(optarg);
        break;
    case 'n':
        searchConfig.leafSize = atoi(optarg);
        std::cerr << " Set leaf node size = " 
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
  
  if ((detectionsFileName == "") || (trackletsFileName == "")){// || (outputFileName == "")) {
    std::cerr << helpString << std::endl;
    return 1;
  }
  
  if(PRINT_TIMING_INFO) {     
    last = std::clock();
  }
  populateDetVectorFromFile(detectionsFileName, allDets);
  populatePairsVectorFromFile(trackletsFileName, allTracklets);
  
  if(PRINT_TIMING_INFO) {     	  
    dif = timeElapsed(last);
    std::cout << "Reading input took " << std::fixed << std::setprecision(10) 
	      <<  dif  << " seconds." <<std::endl;     
  }
  
  
  if(PRINT_TIMING_INFO) {     
    last = std::clock();
  }
  

  /*****************************************************
   * Master node runs linktracklets recursive algorithm
   * and assigns brute force work to slave nodes
   *****************************************************/
  if( rank == 0){
    //run linktracklets program
    std::cout << "Rank " << rank << " calling linkTracklets at " << std::clock() << "." << std::endl;
    linkTracklets(allDets, allTracklets, searchConfig, /*rank,*/ numProcessors);
    std::cout << "Master returned from linkTracklets at " << std::clock() << " and got " << resultTracks.size() << " tracks." << std::endl;
    
    if(PRINT_TIMING_INFO) {     
      dif = timeElapsed (last);
      std::cout << "linking took " << std::fixed << std::setprecision(10) <<  dif 
		<< " seconds."<<std::endl;     
    }
  }
  /*************************************************************
   * Worker nodes wait in a loop to receive tasks for processing
   *************************************************************/
  else if( rank > 0 && rank < numProcessors ){
    std::cerr << "Rank " << rank << " calling waitForTask." << std::endl;
    waitForTask(rank, allDets, allTracklets, searchConfig);
    std::cerr << "Rank " << rank << " returned from waitForTask." << std::endl;
  }
  
  /*
   * Terminate MPI runtime environment
   */
  MPI_Finalize();
  
  return 0;
}

#include <stdlib.h>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <getopt.h>
#include <iomanip>


#include "linkTracklets.h"
#include "../fileUtils.h"



#define PRINT_TIMING_INFO true

#ifdef PRINT_TIMING_INFO
#include <ctime>
#endif

double timeElapsed(clock_t priorEvent)
{
     return ( std::clock() - priorEvent ) / (double)CLOCKS_PER_SEC;
}



void writeResults(std::string outFileName, 
		  const std::vector<Detection> &allDets,
		  const std::vector<Tracklet> &allTracklets,
		  const std::vector<Track> & tracks) 
{
     std::ofstream outFile;
     outFile.open(outFileName.c_str());
     for (unsigned int i = 0; i < tracks.size(); i++) {
	  std::set<unsigned int>::const_iterator detIter;
	  const Track* curTrack = &(tracks.at(i));
	  for (detIter = curTrack->componentDetectionIndices.begin();
	       detIter != curTrack->componentDetectionIndices.end();
	       detIter++) {
	       outFile << *detIter << " ";
	  }
	  outFile << std::endl;
     }
     outFile.close();
}



int main(int argc, char* argv[])
{

     std::string helpString = 
	  "Usage: linkTracklets -d <detections file> -t <tracklets file> -o <output (tracks) file>";

     static const struct option longOpts[] = {
	  { "detectionsFile", required_argument, NULL, 'd' },
	  { "trackletsFile", required_argument, NULL, 't' },
	  { "outputFile", required_argument, NULL, 'o' },
	  { "detectionErrorThresh", required_argument, NULL, 'e'},
	  { "velocityErrorThresh", required_argument, NULL, 'v'},
	  { "maxDecAcceleration", required_argument, NULL, 'D'},
	  { "maxRAAcceleration", required_argument, NULL, 'R'},
	  { "help", no_argument, NULL, 'h' },
	  { NULL, no_argument, NULL, 0 }
     };  
     
     
     std::stringstream ss;
     std::string detectionsFileName = "";
     std::string trackletsFileName = "";
     std::string outputFileName = "";

     linkTrackletsConfig searchConfig; 
     
     int longIndex = -1;
     const char *optString = "d:t:o:e:v:D:R:h";
     int opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
     while( opt != -1 ) {
	  switch( opt ) {
	  case 'd':	       
	       /*ss << optarg; 
		 ss >> detectionsFileName;*/
	       detectionsFileName = optarg;
	       break;
	  case 't':
	       /*ss << optarg;
		 ss >> trackletsFileName; */
	       trackletsFileName = optarg;
	       break;
	  case 'o':
	       /*ss << optarg;
		 ss >> outputFileName; */
	       outputFileName = optarg;
	       break;
	  case 'e':
	       /*ss << optarg;
		 ss >> outputFileName; */
	       searchConfig.detectionLocationErrorThresh = atof(optarg);
	       break;
	  case 'v':
	       /*ss << optarg;
		 ss >> outputFileName; */
	       searchConfig.velocityErrorThresh = atof(optarg);
	       break;
	  case 'D':
	       /*ss << optarg;
		 ss >> outputFileName; */
	       searchConfig.maxDecAccel = atof(optarg);
	       break;
	  case 'R':
	       /*ss << optarg;
		 ss >> outputFileName; */
	       searchConfig.maxRAAccel = atof(optarg);
	       break;
	  case 'h':
	       std::cout << helpString << std::endl;
	       return 0;
	  default:
	       break;
	  }
	  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
     }

     if ((detectionsFileName == "") || (trackletsFileName == "") || (outputFileName == "")) {
	  std::cerr << helpString << std::endl;
	  return 1;
     }

     std::vector<Detection> allDets;
     std::vector<Tracklet> allTracklets;
     std::vector<Track> resultTracks;

     clock_t last;
     double dif;
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

     resultTracks = linkTracklets(allDets, allTracklets, searchConfig);

     if(PRINT_TIMING_INFO) {     
	  dif = timeElapsed (last);
	  std::cout << "linking took " << std::fixed << std::setprecision(10) <<  dif 
		    << " seconds."<<std::endl;     
     }

     if(PRINT_TIMING_INFO) {     
	  last = std::clock();
     }

     writeResults(outputFileName, allDets, allTracklets, resultTracks);

     if(PRINT_TIMING_INFO) {     	  
	  dif = timeElapsed(last);
	  std::cout << "Writing output took " << std::fixed << std::setprecision(10) 
		    <<  dif  << " seconds." <<std::endl;     
     }

     return 0;

    	    
}

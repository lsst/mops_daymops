#include "detectionProximity.h"








void writeResults(std::string, std::vector<std::pair<unsigned int, unsigned int> >);




int main(int argc, char *argv[]){

  std::string dataDetections, queryDetections, outFile;
  double maxDist = 1.0, maxBright = 1.0, maxTime = 1.0;
  outFile = "results.txt";
  
  
  if(argc < 3){
    std::cout << "Usage: detectionProximity -d <data detections> -q "
	      << "<query detections> -o <output file> -t <distance threshold> "
	      << "-b <brightness threshold> -e <time threshold>" << std::endl;
    exit(1);
  }

  static const struct option longOpts[] = {
    { "inFile1", required_argument, NULL, 'd' },
    { "inFile2", required_argument, NULL, 'q' },
    { "outFile", required_argument, NULL, 'o' },
    { "distThresh", required_argument, NULL, 't' },
    { "brightThresh", required_argument, NULL, 'b' },
    { "timeThresh", required_argument, NULL, 'e' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
  };


  int longIndex = -1;
  const char *optString = "d:q:o:t:b:e:h";
  int opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  while( opt != -1 ) {
    switch( opt ) {
    case 'd':
      dataDetections = optarg; 
      break;
    case 'q':
      queryDetections = optarg;
      break;
    case 'o':
      outFile = optarg;
      break;
    case 't':
      maxDist = atof(optarg);
      break;
    case 'b':
      maxBright = atof(optarg);
      break;
    case 'e':
      maxTime = atof(optarg);
      break;
    case 'h':
      std::cout << "Usage: detectionProximity -d <data detections> -q "
	   << "<query detections> -o <output file> -t <distance threshold> "
	   << "-b <brightness threshold> -e <time threshold>" << std::endl;
      exit(0);
    default:
      break;
    }
    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  }


  //open files for reading
  std::ifstream dataDets(dataDetections.c_str());
  std::ifstream queryDets(queryDetections.c_str());

  //used for populating detections vector
  std::vector<MopsDetection> myDataPoints;
  std::vector<MopsDetection> myQueryPoints;
  collapseTracklets::TrackletCollapser myTC; 
  myTC.populateDetVectorFromFile(dataDets, myDataPoints);
  myTC.populateDetVectorFromFile(queryDets, myQueryPoints);

  std::vector<std::pair<unsigned int, unsigned int> > results;
  results = detectionProximity(myQueryPoints,
			       myDataPoints,
			       maxDist,
			       maxBright,
			       maxTime);							

  writeResults(outFile, results);
  
  return(0);
}




/*
 *
 */
void writeResults(std::string outFile, std::vector<std::pair<unsigned int, unsigned int> > results){
  
  std::ofstream writeFile;
  writeFile.open(outFile.c_str());
  
  //write results to file
  for(unsigned int i=0; i<results.size(); i++){
    writeFile << results.at(i).first << " ";
    writeFile << results.at(i).second << std::endl;
  }
  
  writeFile.close();
}


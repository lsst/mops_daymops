// -*- LSST-C++ -*-
/* File: orbitProximity.cc
 * Author: Matthew Cleveland + Jon Myers
 * Purpose: 
 */

#include "lsst/mops/daymops/orbitProximity/orbitProximity.h"


//internal declaration, used only for this file
void writeResults(std::string, std::vector<std::pair<unsigned int, unsigned int> >);

std::vector<lsst::mops::Orbit> populateOrbitVec(std::string);








// main, for calling from command-line
int main(int argc, char *argv[]){

  std::string dataOrbitsFile, queryOrbitsFile, outFile = "results.txt";
  double perihelion = 0.01, eccentricity = .01, inclination = .1; 
  double perihelionArg = 1.0, longitude = 1.0, perihelionTime = .1;
  
  
  if(argc < 3){
      std::cout << "Usage: orbitProximity -d <data orbits> -q "
		<< "<query orbits> -o <output file> -p <perihelion> "
		<< "-e <eccentricity> -i <inclination> -a <perihelion argument> " 
		<< "-l <longitude> -t <perihelion time> -h HELP " << std::endl;
      
      exit(1);
  }

  static const struct option longOpts[] = {
    { "inFile1", required_argument, NULL, 'd' },
    { "inFile2", required_argument, NULL, 'q' },
    { "outFile", required_argument, NULL, 'o' },
    { "periThresh", required_argument, NULL, 'p' },
    { "eccThresh", required_argument, NULL, 'e' },
    { "incThresh", required_argument, NULL, 'i' },
    { "periArgThresh", required_argument, NULL, 'a' },
    { "longThresh", required_argument, NULL, 'l' },
    { "timeThresh", required_argument, NULL, 't' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
  };


  int longIndex = -1;
  const char *optString = "d:q:o:p:e:i:a:l:t:h";
  int opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  while( opt != -1 ) {
    switch( opt ) {
    case 'd':
      dataOrbitsFile = optarg; 
      break;
    case 'q':
      queryOrbitsFile = optarg;
      break;
    case 'o':
      outFile = optarg;
      break;
    case 'p':
      perihelion = atof(optarg);
      std::cerr << "perihelion arg: " << perihelion << std::endl;
      break;
    case 'e':
      eccentricity = atof(optarg);
      std::cerr << "eccentricity arg: " << eccentricity << std::endl;
      break;
    case 'i':
      inclination = atof(optarg);
      std::cerr << "inclination arg: " << inclination << std::endl;
      break;
    case 'a':
      perihelionArg = atof(optarg);
      std::cerr << "perihelionArg arg: " << perihelionArg << std::endl;
      break;
    case 'l':
      longitude = atof(optarg);
      std::cerr << "longitude arg: " << longitude << std::endl;
      break;
    case 't':
      perihelionTime = atof(optarg);
      std::cerr << "perihelionTime arg: " << perihelionTime << std::endl;
      break;
    case 'h':
      std::cout << "Usage: orbitProximity -d <data orbits> -q "
		<< "<query orbits> -o <output file> -p <perihelion> "
		<< "-e <eccentricity> -i <inclination> -a <perihelion argument> " 
		<< "-l <longitude> -t <perihelion time> -h HELP " << std::endl;
      exit(0);
    default:
      break;
    }
    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  }


  //populate Orbit vectors
  std::vector<lsst::mops::Orbit> dataOrbits;
  std::vector<lsst::mops::Orbit> queryOrbits;
  dataOrbits = populateOrbitVec(dataOrbitsFile);
  queryOrbits = populateOrbitVec(queryOrbitsFile);

  std::vector<std::pair <unsigned int, unsigned int> > results = 
      orbitProximity(dataOrbits, queryOrbits,
                     perihelion,
                     eccentricity,
                     inclination,
                     perihelionArg,
                     longitude,
                     perihelionTime);
  
  //write results
  writeResults(outFile, results);
  
  return 0;
}








/**************************************************
 * Read data from 'dataFile' and populate the vector
 * 'dataOrbits' with the Orbit information it contains
 **************************************************/
std::vector<lsst::mops::Orbit> populateOrbitVec(std::string dataFile)
{

    std::vector<lsst::mops::Orbit> dataOrbits;
    std::ifstream myFile(dataFile.c_str());
    std::string line;

    unsigned int lineIndex = 0;

    //assuming file opens properly, iterate through all 
    //lines of file and create orbit object for each
    if(myFile.is_open()){

        while(!myFile.eof()){
            getline(myFile, line);

            //create orbit object from file line
            lsst::mops::Orbit myOrbit;
            myOrbit.populateOrbitFromString(line, lineIndex);
      
            //myOrbit.print();

            //add new Orbit to orbits vector
            dataOrbits.push_back(myOrbit);

            //used for orbit line indexing
            lineIndex++;
        }
    }
    else{
        std::cerr << "Unable to open input file " << dataFile 
                  << "." << std::endl;
        exit(1);
    }
 
    return dataOrbits;
}

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


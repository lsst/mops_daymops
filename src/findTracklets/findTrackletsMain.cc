// -*- LSST-C++ -*-
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <getopt.h>
#include <set>
#include <map>
#include <string>
#include <sstream>
#include <math.h>

#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/fileUtils.h"
#include "lsst/mops/daymops/findTracklets/findTracklets.h"



/*****************************************************************
 *The main program.
 *****************************************************************/
int main(int argc, char* argv[])
{
    //list of all detections
    std::vector <lsst::mops::MopsDetection> myDets; 


    // parse command-line options
    // find name of dets file
    // get maximum velocity (and maybe minimum velocity?) from user
    std::string outFileName;
    std::string inFileName;


    double maxVelocity = 2.0;


    if(argc < 2){
        std::cout << "Usage: findTracklets -i <input file> -o <output file> -v <max velocity> -m <min velocity>" << std::endl;
        exit(1);
    }

    static const struct option longOpts[] = {
        { "inFile", required_argument, NULL, 'i' },
        { "outFile", required_argument, NULL, 'o' },
        { "maxVeloctiy", required_argument, NULL, 'v' },
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
    };

    int longIndex = -1;
    const char *optString = "i:o:v:h";
    int opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
        case 'i':
            inFileName = optarg; 
            break;
        case 'o':
            outFileName = optarg;
            break;
        case 'v':
            maxVelocity = atof(optarg);
            break;
        case 'h':
            std::cout << "Usage: findTracklets -i <input file> -o <output file> -v <max velocity> -m <min velocity>" << std::endl;
            exit(0);
        default:
            break;
        }
        opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    }

    std::ifstream detsFile(inFileName.c_str());
    
    populateDetVectorFromFile(detsFile, myDets);
    
    std::vector<lsst::mops::Tracklet> results;
    
    results = lsst::mops::findTracklets(myDets, maxVelocity);
    
    
    //print results
    std::ofstream outStream;
    outStream.open(outFileName.c_str());

    for(unsigned int i=0; i<results.size(); i++){
        std::set<unsigned int>::const_iterator detIter;
        for (detIter = results.at(i).indices.begin();
             detIter != results.at(i).indices.end();
             detIter++) {
	    outStream << *detIter << " " ;
        }
        outStream << std::endl;
    }
  
    outStream.close();
}


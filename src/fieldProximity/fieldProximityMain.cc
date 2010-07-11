// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
#include <stdlib.h>
#include <algorithm>
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

#include "../KDTree.h"
#include "../TrackletCollapser.h"
#include "Field.h"
#include "fieldProximity.h"


/*****************************************************************
 *The main function of this file.
 *****************************************************************/

int main(int argc, char *args[])
{

    std::string fieldsFile, tracksFile, outFile = "results.txt";
    double maxDist = 1.0;
  
    std::vector<std::string> trackIDs, fieldIDs;
    std::vector<double> fieldTimes;
  
    if(argc < 3){
        std::cout << "Usage: fieldProximity -f <fields file> -t "
                  << "<tracks file> -o <output file> -r <threshold degrees> "
                  << std::endl;
        exit(1);
    }

    static const struct option longOpts[] = {
        { "inFile1", required_argument, NULL, 'f' },
        { "inFile2", required_argument, NULL, 't' },
        { "outFile", required_argument, NULL, 'o' },
        { "distThresh", required_argument, NULL, 'r' },
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
    };


    int longIndex = -1;
    const char *optString = "f:t:o:r:h";
    int opt = getopt_long( argc, args, optString, longOpts, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
        case 'f':
            fieldsFile = optarg; 
            break;
        case 't':
            tracksFile = optarg;
            break;
        case 'o':
            outFile = optarg;
            break;
        case 'r':
            maxDist = atof(optarg);
            break;
        case 'h':
            std::cout << "Usage: detectionProximity -d <data detections> -q "
                      << "<query detections> -o <output file> -t <distance threshold> "
                      << "-b <brightness threshold> -e <time threshold>" << std::endl;
            exit(0);
        default:
            break;
        }
        opt = getopt_long( argc, args, optString, longOpts, &longIndex );
    }
  
  
    //build KDTree from fields file
    //detection vectors, each vector of unique MJD
    std::vector< std::vector<double> > trackSets; 

    //vector of RA and dec pairs for later searching
    std::vector<Field> queryFields; 

    //a vector of Field vectors, each outer vector contains
    //a vector of Fields of a specific time
    std::vector<std::vector<Field> > tracksByTime;

    //link MJD to KDTree of unique MJD detections
    std::map<double, KDTree::KDTree<int> > myTreeMap; 

    //create a vector of all the tracks listed in the tracks
    //file
    std::vector<Field> allFields = populateFieldsVec(fieldsFile, 
                                                     &fieldTimes,
                                                     &fieldIDs);

    std::vector<Field> allTracks = buildTrackVector(tracksFile, 
                                                     &trackIDs);

    //call fieldProximity
    fielProximity(allFields, allTracks);
    

    
    //print results
    printResults(results, fieldIDs, trackIDs, outFile);
    
    exit(EXIT_SUCCESS);
}


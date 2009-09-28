// -*- LSST-C++ -*-
/* jonathan myers */


#include <stdlib.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <getopt.h>
#include <map>


#include "linkTracklets.h"
#include "../fileUtils.h"
#include "../collapseTrackletsAndPostfilters/TrackletCollapser.h"




/*
  this is just a tool for looking through a pair of dets/tracklets
  files and finding objects which meet requirements for being found by
  linkTracklets.
*/


typedef double MJD;

MJD minMJD(std::set<MJD> mjds)
{
    MJD least;
    bool foundOne = false;
    std::set<MJD>::const_iterator mIter;
    for (mIter = mjds.begin(); mIter != mjds.end(); mIter++) {
        if ((!foundOne) || (*mIter < least)) {
            least = *mIter;
            foundOne = true;
        }
    }
    return least;
}


MJD maxMJD(std::set<MJD> mjds)
{
    MJD max;
    bool foundOne = false;
    std::set<MJD>::const_iterator mIter;
    for (mIter = mjds.begin(); mIter != mjds.end(); mIter++) {
        if ((!foundOne) || (*mIter > max)) {
            max = *mIter;
            foundOne = true;
        }
    }
    return max;
}


bool trackletIsCorrect(const std::vector<Detection> & allDets,
                       Tracklet t) 
{
    std::string inferredName = "";
    std::set<unsigned int>::const_iterator detIter;


    for (detIter = t.indices.begin(); detIter != t.indices.end(); detIter++) {

        std::string detObjName = allDets.at(*detIter).getObjName();
        if (inferredName == "") {
            inferredName = detObjName;
        }


        if (detObjName == "NS") { 
            return false; }
        if (detObjName != inferredName) {
            return false;
        }

    }
    return true;
}



MJD firstDetectionTime(const std::vector<Detection> & allDets,
                       const Tracklet &t)
{
    bool foundOne = false;
    MJD minSoFar = 0;
    std::set<unsigned int>::const_iterator detIter;
    for (detIter = t.indices.begin(); detIter != t.indices.end(); detIter++) {
        double curMJD = allDets.at(*detIter).getEpochMJD();
        if ((!foundOne) || (curMJD < minSoFar)) {
            foundOne = true;
            minSoFar = curMJD;
        }
    }
    return minSoFar;
}



void findLinkableObjects(const std::vector<Detection> & allDets, 
                         const std::vector<Tracklet> &allTracklets, 
                         linkTrackletsConfig searchConfig, 
                         std::vector<std::string> &findableObjectNames) 
{

    //consider only detections which come from correct tracklets
    std::map<std::string, std::set<MJD> > nameToObsTimesMap;
    std::map<std::string, std::set<unsigned int> > nameToTrackletsMap;
    
    //get mappings from name to all observation times and name to all component tracklets.
    for (unsigned int trackletI = 0; trackletI != allTracklets.size(); trackletI++) {

        std::set<unsigned int>::const_iterator detIter;
        const Tracklet* curTracklet = &(allTracklets.at(trackletI));
        if (trackletIsCorrect(allDets, *curTracklet)) {
                
            for (detIter = curTracklet->indices.begin(); detIter != curTracklet->indices.end(); detIter++) {
                    
                std::string objName = allDets.at(*detIter).getObjName();
                
                nameToObsTimesMap[objName].insert(allDets.at(*detIter).getEpochMJD());
                nameToTrackletsMap[objName].insert(trackletI);
                
            }
        }
    }


    // iterate through each object and see if its tracklets are legal.
    std::map<std::string, std::set<unsigned int> >::const_iterator objAndTrackletIter;
    for (objAndTrackletIter = nameToTrackletsMap.begin(); 
         objAndTrackletIter != nameToTrackletsMap.end();
         objAndTrackletIter++) {

        std::set<MJD> times = nameToObsTimesMap[objAndTrackletIter->first];

        std::set<MJD> trackletStartTimes;

        std::set<unsigned int>::const_iterator trackletIDIter;
        for (trackletIDIter = objAndTrackletIter->second.begin();
             trackletIDIter != objAndTrackletIter->second.end();
             trackletIDIter++) {
            trackletStartTimes.insert(firstDetectionTime(allDets, allTracklets.at(*trackletIDIter)));
        }
        
        MJD firstTrackletStartTime = minMJD(trackletStartTimes);
        MJD lastTrackletStartTime  = maxMJD(trackletStartTimes);
        
        //check if endpoint separation is sufficient
        if (fabs(lastTrackletStartTime - firstTrackletStartTime) > searchConfig.minEndpointTimeSeparation) {
            
            // collect all detection times and tracklets that could be used
            std::set<MJD> idealTrackDetTimes; 
            std::set<unsigned int> idealTrackTrackletIndices;
            
            // for each tracklet, if it is the first or last tracklet 
            // or a valid support tracklet, add to the idealTrack times and indices.
 
            for (trackletIDIter = objAndTrackletIter->second.begin();
                 trackletIDIter != objAndTrackletIter->second.end();
                 trackletIDIter++) {

                MJD firstObsTime = firstDetectionTime(allDets, allTracklets.at(*trackletIDIter));

                if ((firstObsTime == firstTrackletStartTime) || (firstObsTime == lastTrackletStartTime)) {

                    // this is an ideal endpoint
                    idealTrackTrackletIndices.insert(*trackletIDIter);
                    const Tracklet* curTracklet = &(allTracklets.at(*trackletIDIter));

                    std::set<unsigned int>::const_iterator detIter;
                    for (detIter = curTracklet->indices.begin();
                         detIter != curTracklet->indices.end();
                         detIter++) {
                        idealTrackDetTimes.insert(allDets.at(*detIter).getEpochMJD());
                    }
                }
                else if ((firstObsTime - firstTrackletStartTime > searchConfig.minSupportToEndpointTimeSeparation)
                         && 
                         (lastTrackletStartTime - firstObsTime  > searchConfig.minSupportToEndpointTimeSeparation)) {
                    // this is a valid support tracklet.
                    
                    idealTrackTrackletIndices.insert(*trackletIDIter);
                    const Tracklet* curTracklet = &(allTracklets.at(*trackletIDIter));

                    std::set<unsigned int>::const_iterator detIter;
                    for (detIter = curTracklet->indices.begin();
                         detIter != curTracklet->indices.end();
                         detIter++) {
                        idealTrackDetTimes.insert(allDets.at(*detIter).getEpochMJD());
                    }
                }
            }
           
            // we now have the set of all compatible detections/tracklet. check to see if there are enough.
            if ((idealTrackDetTimes.size() > searchConfig.minDetectionsPerTrack) 
                && 
                (idealTrackTrackletIndices.size() - 2 > searchConfig.minSupportTracklets)){
                // track is valid - add object name to output
                findableObjectNames.push_back(objAndTrackletIter->first);
            }
        }
    }
        
}





int main(int argc, char** argv) 
{
     std::string helpString = 
	  "Usage: linkTracklets -d <detections file> -t <tracklets file> -o <output (tracks) file>";
     
     static const struct option longOpts[] = {
	  { "detectionsFile", required_argument, NULL, 'd' },
	  { "trackletsFile", required_argument, NULL, 't' },
	  { "outputFile", required_argument, NULL, 'o' },
	  { "help", no_argument, NULL, 'h' },
	  { NULL, no_argument, NULL, 0 }
     };  
     
     
     std::stringstream ss;
     std::string detectionsFileName = "";
     std::string trackletsFileName = "";
     std::string outputFileName = "";
     
     int longIndex = -1;
     const char *optString = "d:t:o:h";
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


     populateDetVectorFromFile(detectionsFileName, allDets);
     populatePairsVectorFromFile(trackletsFileName, allTracklets);
     
     linkTrackletsConfig searchConfig; 

     std::vector<std::string> findableObjectNames;

     findLinkableObjects(allDets, allTracklets, searchConfig, findableObjectNames);

     // write to disk
     std::ofstream outFile;
     outFile.open(outputFileName.c_str());
     for (unsigned int i = 0; i < findableObjectNames.size(); i++) {
         outFile << findableObjectNames.at(i) << std::endl;
     }
     outFile.close();

}

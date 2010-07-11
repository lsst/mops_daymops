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
 

/*
 * jmyers 9/20/09
 * 
 */

#include <istream>
#include <sstream>

#include "lsst/mops/fileUtils.h"


namespace lsst {
    namespace mops {






bool isSane(unsigned int detsSize, const std::vector<Tracklet> *pairs) {
    std::vector<Tracklet>::const_iterator tIter;
    std::set<unsigned int>::const_iterator iIter;
    for (tIter = pairs->begin(); tIter != pairs->end(); tIter++) {
        for (iIter = tIter->indices.begin(); iIter != tIter->indices.end(); 
             iIter++) {
            if (*iIter >= detsSize) {
                return false;
            }
        }
    }
    return true;
}








void writeTrackletsToOutFile(const std::vector<Tracklet> * tracklets, std::ofstream &outFile){
     std::vector<Tracklet>::const_iterator trackletIter;
     std::set<unsigned int>::iterator indicesIter;
     for (trackletIter = tracklets->begin(); trackletIter != tracklets->end(); trackletIter++) {
	  for (indicesIter = trackletIter->indices.begin(); 
	       indicesIter != trackletIter->indices.end();
	       indicesIter++) {
	       outFile << *indicesIter;
	       outFile << " ";
	  }
	  outFile << std::endl;
     }
}







void populateDetVectorFromFile(std::ifstream &detsFile, std::vector <MopsDetection> &myDets) {
     std::string line;
     MopsDetection tmpDet;
     line.clear();
     std::getline(detsFile, line);
     while (detsFile.fail() == false) {
	  tmpDet.fromMITIString(line);
	  myDets.push_back(tmpDet);
	  line.clear();
	  std::getline(detsFile, line);
     }
}






void populatePairsVectorFromFile(std::ifstream &pairsFile, 
						    std::vector<Tracklet> &pairsVector) {

     Tracklet tmpPair;
     std::string line;
     int tmpInt = -1;
     tmpPair.indices.clear();
     tmpPair.isCollapsed = false;
     line.clear();
     std::getline(pairsFile, line);
     while (pairsFile.fail() == false) {
	  std::istringstream ss(line);
	  ss.exceptions(std::ifstream::failbit | std::ifstream::badbit);    
	  while (!ss.eof()) {
	       try {
                    ss >> tmpInt;
                    ss >> std::ws;
	       }
	       catch (...) {
                    throw LSST_EXCEPT(InputFileFormatErrorException, "Improperly-formatted pairs file.\n");
	       }
	       tmpPair.indices.insert(tmpInt);
	  }
	  if (tmpPair.indices.size() < 2) {
	       throw LSST_EXCEPT(InputFileFormatErrorException, "EE: CollapseTracklets: pairs in pairs file must be length >= 2!\n");
	  }
	  pairsVector.push_back(tmpPair);
	  tmpPair.indices.clear();
	  tmpPair.isCollapsed = false;
	  line.clear();
	  std::getline(pairsFile, line);   
     }
}
  


/* these overloaded versions are mainly for SWIG, since Python file objects != fstreams */
void writeTrackletsToOutFile(const std::vector<Tracklet> * tracklets, std::string outFileName)
{
     std::ofstream outFile;
     outFile.open(outFileName.c_str());
     if (!outFile.is_open()) {
	  throw LSST_EXCEPT(FileException,
			    "Failed to open output file " + outFileName + " - do you have permission?\n");
     }
     writeTrackletsToOutFile(tracklets, outFile);
}
    
void populateDetVectorFromFile(std::string detsFileName, std::vector <MopsDetection> &myDets)
{
     std::ifstream detsFile;
     detsFile.open(detsFileName.c_str());
     if (!detsFile.is_open()) {
	  throw LSST_EXCEPT(FileException,
			    "Failed to open dets file " + detsFileName + " - does this file exist?\n");
     }
     populateDetVectorFromFile(detsFile, myDets);
        
}
    
void populatePairsVectorFromFile(std::string pairsFileName,
						    std::vector <Tracklet> &pairsVector)
{
     std::ifstream pairsFile;
     pairsFile.open(pairsFileName.c_str());
     if (!pairsFile.is_open()) {
	  throw LSST_EXCEPT(FileException,
			    "Failed to open pairs file " + pairsFileName + " - does this file exist?\n");
     }
     populatePairsVectorFromFile(pairsFile, pairsVector);
        
}

    }} // close lsst::mops

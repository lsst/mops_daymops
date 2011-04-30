// -*- LSST-C++ -*-
/* jonathan myers 
   4/08/10
*/

//#include <string>
#include <iomanip>
#include <ctime>
#include <time.h>

#include "lsst/mops/Track.h"
#include "lsst/mops/Exceptions.h"
#include "lsst/mops/TrackVector.h"


namespace lsst { 
namespace mops {


void TrackVector::populateFromFile(std::string fileName, 
                                   const std::vector<MopsDetection> &allDets)
{
    std::ifstream trackFile;
    trackFile.open(fileName.c_str());
    Track *tmpTrack;
    std::string line;
    unsigned int tmpInt = 0;
    std::getline(trackFile, line);
    while (trackFile.fail() == false) {
        tmpTrack = new Track;
        std::istringstream ss(line);
        ss.exceptions(std::ifstream::failbit | std::ifstream::badbit);    
        while (!ss.eof()) {
            try {
                ss >> tmpInt;
                ss >> std::ws;
            }
            catch (...) {
                delete tmpTrack;
                throw LSST_EXCEPT(InputFileFormatErrorException, 
                                  "Improperly-formatted track file.\n");
            }
            if (tmpInt >= allDets.size()) {
                delete tmpTrack;
                throw LSST_EXCEPT(InputFileFormatErrorException,
                                  "Got index larger than size of detections vector. Could this file be in the wrong format (e.g. DiaSourceIDs instead of indices) or is this the wrong detections file?");
            }
            tmpTrack->addDetection(tmpInt, allDets);
        }
        contents.push_back(*tmpTrack);
        line.clear();
        std::getline(trackFile, line);   
        delete tmpTrack;
    }
}






void TrackVector::writeSet(std::ofstream *outFile, const std::set <unsigned int> &s)
{
    std::set<unsigned int>::const_iterator outIter;
    for (outIter = s.begin(); outIter != s.end(); outIter++) {
        *outFile << *outIter  << " ";
    }
}







void TrackVector::writeTracksToFile(std::string outFileName)
{

     std::ofstream outFile;
     outFile.open(outFileName.c_str());
     for (unsigned int i = 0; i < contents.size(); i++) {
         writeSet(&outFile, this->at(i)->getComponentDetectionDiaIds());
         outFile << "\n";
     }
}






void TrackVector::writeTracksAndStatsToFile(std::string outFileName,
                                            const std::vector<MopsDetection> &allDets) 
{
     std::ofstream outFile;
     outFile.open(outFileName.c_str());
     for (unsigned int i = 0; i < contents.size(); i++) {
         Track* curTrack = this->at(i);
         writeSet(&outFile, curTrack->getComponentDetectionDiaIds());
         curTrack->calculateBestFitQuadratic(allDets, true);
         double epoch, ra0, raV, raAcc, dec0, decV, decAcc;
         curTrack->getBestFitQuadratic(epoch, ra0, raV, raAcc, 
                                      dec0, decV, decAcc);
         outFile << " epoch:" << std::scientific << epoch
                 << " Best-fit RA p0, v, acc: " 
                 << ra0 << " " 
                 << raV << " " 
                 << raAcc << " " 
                 << " Best-fit Dec p0, v, acc: " 
                 << dec0 << " " 
                 << decV << " " 
                 << decAcc  
                 << " Underlying object: " 
                 << curTrack->getObjectId(allDets) 
                 << " Chi squared prob: Ra: " 
                 << curTrack->getProbChisqRa() 
                 << " Dec: "
                 << curTrack->getProbChisqDec()
                 << "\n";
     }
}




}} // close namespace lsst::mops

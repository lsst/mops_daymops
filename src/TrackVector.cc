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


}} // close namespace lsst::mops

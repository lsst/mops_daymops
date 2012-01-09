// -*- LSST-C++ -*-

#include "lsst/mops/Track.h"
#include "lsst/mops/TrackVector.h"
#include "lsst/mops/fileUtils.h"
#include "lsst/mops/daymops/linkTracklets/linkTracklets.h"

#include <string>
#include <iostream>


namespace lsst { namespace mops {


void writeSet(std::ofstream *outFile, const std::set <unsigned int> &s)
{
    std::set<unsigned int>::const_iterator outIter;
    for (outIter = s.begin(); outIter != s.end(); outIter++) {
        *outFile << *outIter  << " ";
    }
}


void processTracks(std::string fileName, std::string outFileName,
                   const std::vector<MopsDetection> &allDets)
{
    // open infile
    std::ifstream trackFile;
    trackFile.open(fileName.c_str());
    Track *tmpTrack = NULL;
    std::string line;
    unsigned int tmpInt = 0;
    std::getline(trackFile, line);


    // open outfile
    std::ofstream outFile;
    outFile.open(outFileName.c_str());

    // try to read a track.
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

        if (tmpTrack != NULL) {
            // got a track. calc stats and write to file.
            Track curTrack = *tmpTrack;
            writeSet(&outFile, curTrack.getComponentDetectionDiaIds());
            curTrack.calculateBestFitQuadratic(allDets, true);
            double epoch, ra0, raV, raAcc, dec0, decV, decAcc;
            curTrack.getBestFitQuadratic(epoch, ra0, raV, raAcc, 
                                         dec0, decV, decAcc);
            double chisqR = curTrack.getProbChisqRa();
            double chisqD = curTrack.getProbChisqDec();
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
                    << curTrack.getObjectId(allDets) << " "
                    << "Chi Squared prob (Ra, Dec, Ra*Dec) : "
                    << chisqR << " "
                    << chisqD << " "
                    << chisqR*chisqD 
                 << "\n";

            delete tmpTrack;
            tmpTrack = NULL;
            
        }
        
        
        
        line.clear();
        std::getline(trackFile, line);   
    }


}

int doIt(int argc, char** argv) {


     if (argc != 4) {
	  std::string usage("USAGE: readTracksWriteStats <inDets> <inTracks_byIndices> <outFile> \n");
	  std::cout << usage;
	  std::cerr << usage;
	  return 1;
     }

     lsst::mops::linkTrackletsConfig searchConfig;

     std::string inDetsFile(argv[1]);
     std::string inFile(argv[2]);
     std::string outFile(argv[3]);
     
     TrackVector myTv;
     std::vector<MopsDetection> allDets;
     std::cout << "Reading dets...\n";
     populateDetVectorFromFile(inDetsFile, allDets, searchConfig.defaultAstromErr);
     for (uint i = 0; i < allDets.size(); i++) {
         allDets[i].setObservatoryLocation(searchConfig.obsLat, searchConfig.obsLong);
         allDets[i].calculateTopoCorr();
     }
     std::cout << "Reading input tracks (as indices) and writing output...\n";

     processTracks(inFile, outFile, allDets);

     std::cout << "Done.\n";

     return 0;
}

}} // close lsst::mops




int main(int argc, char** argv) 
{
     return lsst::mops::doIt(argc, argv);

}

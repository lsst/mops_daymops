// -*- LSST-C++ -*-


/* jmyers 8/18/08 
 */

#include <iomanip>
#include <sstream>

#include <unistd.h>
#include <getopt.h>

#include "lsst/mops/fileUtils.h"
#include "lsst/mops/rmsLineFit.h"

namespace lsst {
    namespace mops {    







    /* purifyTracklet:  assuming *t is an allocated tracklet, and allDets is an allocated vector of
     * MopsDetections into which t's indices are indeed indices.
     * 
     * return a corresponding tracklet, possibly empty, s.t. the tracklet's max
     * RMS is less than average magnitude * maxRMSm + maxRMSb.  This tracklet
     * contains a subset of the detections associated with *t.  If t actually
     * contains multiple tracklets, this will not help; but it will probably
     * help you find one tracklet, anyway.
     */
    int rmsPurifyMain(int argc, char** argv) {
        clock_t start = std::clock();

        std::string USAGE("USAGE: purifyTracklets --detsFile <detections file> --pairsFile <tracklets (pairs) file) --maxRMS --outFile <output tracklets (pairs) file>");
        char* pairsFileName = NULL;
        char* detsFileName = NULL;
        char* outFileName = NULL;

        std::vector <Tracklet> trackletsVector;
        std::vector <MopsDetection> detsVector;
        
        double maxRMS = .001;
        unsigned int minObs = 2;

        static const struct option longOpts[] = {
            { "pairsFile", required_argument, NULL, 'p' },
            { "detsFile", required_argument, NULL, 'd' },
            { "outFile", required_argument, NULL, 'o' },
            { "maxRMS", required_argument, NULL, 'm'},
            { "minobs", required_argument, NULL, 'n'},
            { "help", no_argument, NULL, 'h' },
            { NULL, no_argument, NULL, 0 }
        };

        int longIndex = -1;
        const char* optString = "p:d:o:m:h";        
        int opt = getopt_long( argc, argv, optString, longOpts, &longIndex );        
        std::stringstream ss;

        while( opt != -1 ) {
            switch( opt ) {
            case 'p':
                pairsFileName = optarg; 
                break;
                
            case 'd':
                detsFileName = optarg ; 
                break;
                
            case 'o':
                outFileName = optarg;
                break;
                
            case 'm':
                ss.clear();
                ss << optarg;
                ss >> maxRMS;
                break;

            case 'n':
                ss.clear();
                ss << optarg;
                ss >> minObs;
                break;
                
            case 'h':   /* fall-through is intentional */
            case '?':
                std::cout << "got request for help " << std::endl;
                std::cout<<USAGE<<std::endl;
                return 0;
                break;
            default:
                throw LSST_EXCEPT(ProgrammerErrorException,
                                  "EE: Unexpected programmer error in options parsing\n");
                break;
            }        
            opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
        }

        if ((pairsFileName == NULL) || (detsFileName == NULL) || (outFileName == NULL)) {
            throw LSST_EXCEPT(CommandlineParseErrorException,
                              "Did not get required parameters \n\n"
                              + USAGE + "\n");
        }
        std::ifstream pairsFile(pairsFileName);
        std::ifstream detsFile(detsFileName);
        std::ofstream outFile(outFileName);

        if (!detsFile.is_open()) {
            throw LSST_EXCEPT(FileException,
                              "Failed to open dets file " + std::string(detsFileName) + " - does this file exist?\n");                          
        }
        if (!pairsFile.is_open()) {
            throw LSST_EXCEPT(FileException,
                              "Failed to open pairs file " + std::string(pairsFileName) + " - does this file exist?\n");                          
        }
        if (!outFile.is_open()){
            throw LSST_EXCEPT(FileException,
                              "Failed to open output file " + 
                              std::string(outFileName) + " - do you have write permissions?\n");
        }

        std::cout << "purifyTracklets: " << std::endl;
        std::cout << "==========================================" << std::endl;
        std::cout << "Detections file:        " << detsFileName << std::endl;
        std::cout << "Pairs (tracklets) file: " << pairsFileName << std::endl;
        std::cout << "Output file:            " << outFileName << std::endl;
        std::cout << "max RMS:                " <<  maxRMS << std::endl;

        std::cout << "Reading detections file...." << std::endl;
        populateDetVectorFromFile(detsFile, detsVector);
        std::cout << "Reading tracklets (pairs) file...." << std::endl;
        populatePairsVectorFromFile(pairsFile, trackletsVector);
        std::cout << "Done!" << std::endl;
        double dif = lsst::mops::timeElapsed(start);
        std::cout << "Reading input took " << std::fixed << std::setprecision(10) 
                  <<  dif  << " seconds." <<std::endl;             

        if (!isSane(detsVector.size(), &trackletsVector)) {
            throw LSST_EXCEPT(InputFileFormatErrorException, 
                              "EE: Pairs file does not seem to correspond with detections file.\n");
        }

        std::vector<Tracklet> postFilteredTracklets;        
        std::cout << "Doing the filtering..." << std::endl;
        TrackletPurifier myTP;
        
        myTP.purifyTracklets(&trackletsVector, &detsVector, maxRMS, minObs, postFilteredTracklets);
        
        std::cout << "Done. Writing output." << std::endl;
        writeTrackletsToOutFile(&postFilteredTracklets, outFile);

        std::cout << "Done!" << std::endl;
        std::cout << "Completed after " << std::fixed << std::setprecision(10) 
                  <<  dif  << " seconds." <<std::endl;
        printMemUse();
        return 0;
    }






}} // close lsst::mops

int main(int argc, char** argv) {
    return lsst::mops::rmsPurifyMain(argc, argv);
}



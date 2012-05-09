// -*- LSST-C++ -*-
/* jonathan myers */

#include <iomanip>
#include <sstream>

#include <unistd.h>
#include <getopt.h>

#include "lsst/mops/common.h"
#include "lsst/mops/daymops/collapseTrackletsAndPostfilters/collapseTracklets.h"
#include "lsst/mops/fileUtils.h"

namespace lsst {
    namespace mops {


    bool detectionsHaveGT2UniqueTimes(const std::vector<MopsDetection> *detections) {
        std::vector<double> observedMJDs;
        double curMJD = -1;
        std::vector<double>::iterator fIter;
        for (unsigned int i = 0; i < detections->size(); i++) {

            curMJD = (*detections)[i].getEpochMJD();
            bool isInList = false;

            for (fIter = observedMJDs.begin(); 
                 (fIter != observedMJDs.end()) && (isInList == false);
                 fIter++) {
                if (areEqual(*fIter, curMJD)) {
                    isInList = true;
                }
            }

            if (isInList == false) {
                observedMJDs.push_back(curMJD);
            }
            if (observedMJDs.size() > 2) {
                return true;
            }
        }
        /* we never got observedMJDs > 2, return false */
        return false;        
    }


    int collapseTrackletsMain(int argc, char** argv)
    { 
      time_t start = time(NULL);
        std::string USAGE("USAGE: collapseTracklets [options] <detsFile> <pairsFile> <RA Tolerance> <Dec Tolerance> <angular tolerance> <velocity tolerance> <outFile>\n\
options:  \n\
=================================================\n\
--method=<greedy|minimumRMS|bestFit>    \n\
-------------------------------------------------\n\
if greedy, then we choose as many compatible tracklets as possible, as returned by the tree search.  if minimumRMS, we take the results of the tree search and repeatedly find the tracklet which would have the lowest resulting RMS value if added, then add it. If bestFit, we repeatedly choose the tracklet which is closest to the current approximate line first, rather than re-calculating best-fit line for each possible tracklet. (NB: Greedy is default.)\n\
\n\
--useRMSFilt=<true/false>  \n\
--------------------------------------------------\n\
enforce a maximum RMS distance for any tracklet which is the product of collapsing.  default is false.\n\
\n\
--maxRMS=<double>,\n\
--------------------------------------------------\n\
Only used if useRMSfilt == true.  Describes the function for RMS filtering.  Tracklets will not be collapsed unless the resulting tracklet would have RMS <= maxRMSm * average magnitude + maxRMSb.   Defaults are 0. and .001.\n\
");

        std::ifstream detsFile;
        std::ifstream pairsFile;
        std::vector <double> tolerances;
        std::ofstream outFile;
        std::string line;
        std::vector <MopsDetection> detections;
        /* pairs holds indices into detections. */
        std::vector <Tracklet> pairs;
        std::vector<Tracklet> collapsedPairs;        

        /* read parameters from command-line */

        /* read long options, then get required options */

        if (argc < 8) {
            throw LSST_EXCEPT(CommandlineParseErrorException, USAGE + "\n");            
        }

        
        bool useMinimumRMS = false;
	bool useGreedy = true;
        bool useBestFit = false;
        bool useRMSFilt = false;
        double maxRMS = .001;
        
        static const struct option longOpts[] = {
            { "method", required_argument, NULL, 'e' },
            { "useRMSFilt", required_argument, NULL, 'u' },
            { "maxRMS", required_argument, NULL, 'm' },
            { "help", no_argument, NULL, 'h' },
            { NULL, no_argument, NULL, 0 }
        };

        const char* optString = "e:u:m:h:v";
        int longIndex = -1;
        int opt = getopt_long( argc, argv, optString, longOpts, &longIndex );        
        std::stringstream ss;
	bool beVerbose = false;
        while( opt != -1 ) {
            switch( opt ) {
	    case 'v':
                beVerbose = true;
                break;
            case 'e':
                if (std::string(optarg) == std::string("greedy")) {
                    useMinimumRMS = false;                  
		    useBestFit = false;
                }
                else if (std::string(optarg) == std::string("minimumRMS")) {
                    useMinimumRMS = true;
		    useGreedy = false;
                }
                else if (std::string(optarg) == std::string("bestFit")) {
                    useBestFit = true;
		    useMinimumRMS = false;
		    useGreedy = false;
                }
                else {
                    throw LSST_EXCEPT(CommandlineParseErrorException, 
                                      "ERROR: options for method are: greedy, minimumRMS\n\n" +
                                      USAGE + "\n"); 
                }
                break;

            case 'u':
                useRMSFilt = guessBoolFromStringOrGiveErr(optarg, USAGE);
                break;
                
            case 'm':
                ss.clear();
                ss << optarg;
                ss >> maxRMS;
                break;
                
            case 'h':   /* fall-through is intentional */
            case '?':
                std::cout << "got request for help " << std::endl;
                std::cout<<USAGE<<std::endl;
                return 0;
                break;
            default:
                throw LSST_EXCEPT(ProgrammerErrorException, 
                                  "Unexpected programmer error in options parsing\n");
                break;
            }        
            opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
        }

        std::string detsFileName;
        std::string pairsFileName;
        std::string outFileName;

        /* note that by opening these fileswq these here we do allow some
         * redundant code, since we have functions which take a string and open
         * the file by the same name, but we want to check that the parameters
         * are good *before* we do a bunch of work.
         */

        detsFileName = argv[optind];
        detsFile.open(detsFileName.c_str());
        if (!detsFile.is_open()) {
            throw LSST_EXCEPT(FileException,
                              "Failed to open dets file " + detsFileName + " - does this file exist?\n"); 
        }
        optind++;
        pairsFileName = argv[optind];
        pairsFile.open(pairsFileName.c_str());
        if (!pairsFile.is_open()) {
            throw LSST_EXCEPT(FileException,
                              "Failed to open pairs file " + pairsFileName + " - does this file exist?\n");
        }
        optind++;
        for (int i = optind; i < optind + 4; i++) {
            double tmp = 0.;
            std::stringstream ss(argv[i]);
            ss >> tmp;
            tolerances.push_back(tmp);
        }
        optind += 4;
        outFileName = argv[optind];
        outFile.open(outFileName.c_str());
        if (!outFile.is_open()){
            throw LSST_EXCEPT(FileException,
                              "Failed to open output file " + 
                              outFileName + " - do you have write permissions?\n");
        }


        std::cout << "collapseTracklets: " << std::endl << "--------------------------------" << std::endl;
        std::cout << "Method:            ";
        if (useMinimumRMS) {
            std::cout << "minimum RMS";
        }
        if (useGreedy) { 
            std::cout << "greedy" ;
        }
        if (useBestFit) {
            std::cout << "best fit";
        }
        std::cout << std::endl;
        if (useRMSFilt == true) {
            std::cout << "enforce min. RMS:  true" << std::endl;
            std::cout << "maximum RMS:       " << maxRMS << std::endl;
        }
        std::cout << "Dets File:         " << detsFileName << std::endl;
        std::cout << "Pairs File         " << pairsFileName << std::endl;
        std::cout << "Output File:       " << outFileName << std::endl;
        std::cout << "Tolerances: " << std::endl;
        std::cout << "RA0: " << tolerances[0] << ", Dec0: " << tolerances[1] << ", angle: " << 
            tolerances[2] << ", velocity: " << tolerances[3] << std::endl;


        /* read data from files, do collapsing if data supports it, write
         * output */

        populateDetVectorFromFile(detsFile, detections);
        populatePairsVectorFromFile(pairsFile, pairs);
        
        double dif = lsst::mops::timeElapsed(start);
        std::cout << "Reading input took " << std::fixed << std::setprecision(10) 
                  <<  dif  << " seconds." <<std::endl;             

        // TBD: use proper LSST logs, too...
        std::cout << "Found " << detections.size() << " detections\n";
        std::cout << "Found " << pairs.size() << " pairs.\n";
          
        if (!isSane(detections.size(), &pairs)) {
            throw LSST_EXCEPT(InputFileFormatErrorException, 
                              "EE: Pairs file does not seem to correspond with detections file.\n");
        }
          

        if ((detectionsHaveGT2UniqueTimes(&detections)) && (pairs.size() > 1)) {
            std::cout << "Doing the collapsing..." << std::endl;
            doCollapsingPopulateOutputVector(&detections, pairs, tolerances, collapsedPairs, 
                                             useMinimumRMS, useBestFit, useRMSFilt, maxRMS, 
                                             beVerbose);
            std::cout << std::endl << "Done!" << std::endl;

        }
        else {
            /* no need to do any collapsing, there isn't enough data! */
            collapsedPairs = pairs;
        }
        std::cout << "Writing the output to file..." << std::endl;
        /* write the output */
        /*collapsedPairs = removeSubsets(&collapsedPairs);  for now, do this in a separate program.*/
        writeTrackletsToOutFile(&collapsedPairs, outFile);
        /* close all files. */
        detsFile.close();
        pairsFile.close();
        outFile.close();	

        std::cout << "Done!" << std::endl;
        std::cout << "Completed after " << std::fixed << std::setprecision(10) 
                  <<  dif  << " seconds." <<std::endl;
        printMemUse();
        return 0;
    }
  









}} // close lsst::mops



int main(int argc, char** argv) {
        
    return lsst::mops::collapseTrackletsMain(argc, argv);
    
}


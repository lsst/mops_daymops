// -*- LSST-C++ -*-
/* jonathan myers */

#include <iomanip>

#include <unistd.h>
#include <getopt.h>

#include "lsst/mops/fileUtils.h"
#include "lsst/mops/removeSubsets.h"



namespace lsst {
namespace mops {

  extern std::string curTime();

    int removeSubsetsMain(int argc, char** argv) {
      time_t start = time(NULL);

        /* read pairs from a file, do the removal, and write the output */
        std::string USAGE("USAGE: removeSubsets --inFile <input pairfile> --outFile <output pairfile> [--removeSubsets <TRUE/FALSE> --keepOnlyLongest <TRUE/FALSE>]");
        std::ifstream inFile;
        std::ofstream outFile;
        std::vector <Tracklet>* pairsVector = new std::vector<Tracklet>;
        std::vector <Tracklet>* outputVector = new std::vector<Tracklet>;

        char* inFileName = NULL;
        char* outFileName = NULL;
        bool removeSubsets = true;
        bool keepOnlyLongestPerDet = false;
        bool shortCircuit = true;
        bool sortBeforeIntersect = false;
        
        static const struct option longOpts[] = {
            { "inFile", required_argument, NULL, 'i' },
            { "outFile", required_argument, NULL, 'o' },
            { "removeSubsets", required_argument, NULL, 'r' },
            { "keepOnlyLongest", required_argument, NULL, 'k'},
            { "shortCircuit", required_argument, NULL, 'c'},
            { "sortBeforeIntersect", required_argument, NULL, 's'},
            { "help", no_argument, NULL, 'h' },
            { NULL, no_argument, NULL, 0 }
        };

        int longIndex = -1;
        const char* optString = "i:o:r:k:c:s:h";
        int opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
        while( opt != -1 ) {
            switch( opt ) {
            case 'i':
                inFileName = optarg; 
                break;
                
            case 'o':
                outFileName = optarg ; 
                break;
                
            case 'r':
                removeSubsets = guessBoolFromStringOrGiveErr(optarg, USAGE);
                break;
                
            case 'k':
                keepOnlyLongestPerDet = guessBoolFromStringOrGiveErr(optarg, USAGE);
                break;
            case 'c':
                shortCircuit = guessBoolFromStringOrGiveErr(optarg, USAGE);
                break;
            case 's':
                sortBeforeIntersect = guessBoolFromStringOrGiveErr(optarg, USAGE);
                break;
                
            case 'h':   /* fall-through is intentional */
            case '?':
                std::cout<<USAGE<<std::endl;
                return 0;
                break;
            default:
                throw LSST_EXCEPT(ProgrammerErrorException, "Programmer error in parsing of command-line options\n");
                break;
            }        
            opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
        }

        if ((inFileName == NULL) || (outFileName == NULL)) {
            throw LSST_EXCEPT(CommandlineParseErrorException, 
                              "Please specify input and output filenames.\n\n" +
                              USAGE + "\n");
        }


        std::cout << "RemoveSubsets:" << std::endl;
        std::cout << "---------------" << std::endl;
        std::cout << "Input file:                                  " << inFileName << std::endl;
        std::cout << "Output file:                                 " << outFileName << std::endl;
        std::cout << "RemoveSubsets:                               "<< boolToString(removeSubsets)
                  << std::endl;
        std::cout << "Short-circuit if possible:                   "<< boolToString(shortCircuit)
                  << std::endl;
        std::cout << "Sort sets by size before set intersection:   "<< boolToString(sortBeforeIntersect)
                  << std::endl;
        std::cout << "Keep only longest tracklet(s) per detection: "<< 
            boolToString(keepOnlyLongestPerDet) << std::endl;

        inFile.open(inFileName);
        outFile.open(outFileName);
	std::cout << "Reading infile, starting at " << curTime() << std::endl;
        populatePairsVectorFromFile(inFile, *pairsVector);
	std::cout << "Finished reading infile at " << curTime() << std::endl;
        double dif = lsst::mops::timeElapsed(start);
        std::cout << "Reading input took " << std::fixed << std::setprecision(10) 
                  <<  dif  << " seconds." <<std::endl;             

        /* do the actual work */

        if (keepOnlyLongestPerDet == true) {
            std::vector <Tracklet>* tmpVector = 
                new std::vector<Tracklet>;
            putLongestPerDetInOutputVector(pairsVector, *tmpVector); 
            delete pairsVector;
            pairsVector = tmpVector;
        }

        SubsetRemover mySR;
        
        if (removeSubsets == true) {
            mySR.removeSubsetsPopulateOutputVector(pairsVector, 
                                                   *outputVector,
                                                   shortCircuit, 
                                                   sortBeforeIntersect);
        }
        else {
            delete outputVector;
            outputVector = pairsVector;
        }


	std::cout << "Finished filtering, writing output starting at " << curTime() << std::endl;
        writeTrackletsToOutFile(outputVector, outFile);
	outFile.close();
        if (outputVector != pairsVector) {
            delete outputVector;
        }
        delete pairsVector;

        std::cout << "Completed after " << std::fixed << std::setprecision(10) 
                  <<  dif  << " seconds." <<std::endl;

	if (outFile.good()) {
	  std::cout << "Finished successfully finished at " << curTime() << std::endl;
          printMemUse();
	  return 0;
	}
	else {
	  std::cout << "ERROR writing/closing file. Finished at " << curTime() << std::endl;
          printMemUse();
	  return 1;
	}
    }



}} // close lsst::mops

int main(int argc, char** argv) {
    return lsst::mops::removeSubsetsMain(argc, argv);
}

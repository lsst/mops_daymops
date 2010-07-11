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
 
/* jonathan myers */

#include <unistd.h>
#include <getopt.h>

#include "SubsetRemover.h"


namespace ctExcept = collapseTracklets::exceptions;
namespace removeSubsets {

    int removeSubsetsMain(int argc, char** argv) {

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
        
        static const struct option longOpts[] = {
            { "inFile", required_argument, NULL, 'i' },
            { "outFile", required_argument, NULL, 'o' },
            { "removeSubsets", required_argument, NULL, 'r' },
            { "keepOnlyLongest", required_argument, NULL, 'k'},
            { "help", no_argument, NULL, 'h' },
            { NULL, no_argument, NULL, 0 }
        };

        int longIndex = -1;
        const char* optString = "i:o:r:k:h";
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
                removeSubsets = KDTree::Common::guessBoolFromStringOrGiveErr(optarg, USAGE);
                break;
                
            case 'k':
                keepOnlyLongestPerDet = KDTree::Common::guessBoolFromStringOrGiveErr(optarg, USAGE);
                break;
                
            case 'h':   /* fall-through is intentional */
            case '?':
                std::cout<<USAGE<<std::endl;
                return 0;
                break;
            default:
                throw LSST_EXCEPT(ctExcept::ProgrammerErrorException, "Programmer error in parsing of command-line options\n");
                break;
            }        
            opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
        }

        if ((inFileName == NULL) || (outFileName == NULL)) {
            throw LSST_EXCEPT(ctExcept::CommandlineParseErrorException, 
                              "Please specify input and output filenames.\n\n" +
                              USAGE + "\n");
        }


        std::cout << "RemoveSubsets:" << std::endl;
        std::cout << "---------------" << std::endl;
        std::cout << "Input file:                                  " << inFileName << std::endl;
        std::cout << "Output file:                                 " << outFileName << std::endl;
        std::cout << "RemoveSubsets:                               "<< KDTree::Common::boolToString(removeSubsets)
                  << std::endl;
        std::cout << "Keep only longest tracklet(s) per detection: "<< 
            KDTree::Common::boolToString(keepOnlyLongestPerDet) << std::endl;

        collapseTracklets::TrackletCollapser myTC;

        inFile.open(inFileName);
        outFile.open(outFileName);
        myTC.populatePairsVectorFromFile(inFile, *pairsVector);

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
            mySR.removeSubsetsPopulateOutputVector(pairsVector, *outputVector);
        }
        else {
            delete outputVector;
            outputVector = pairsVector;
        }

        myTC.writeTrackletsToOutFile(outputVector, outFile);
        if (outputVector != pairsVector) {
            delete outputVector;
        }
        delete pairsVector;
        return 0;
    }

}

int main(int argc, char** argv) {
    CALL_AND_CATCH_EXCEPTIONS(return removeSubsets::removeSubsetsMain(argc, argv));
}

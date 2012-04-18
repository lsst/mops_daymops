// -*- LSST-C++ -*-
/* jonathan myers */

#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <sstream>

#include "lsst/mops/common.h"
#include "lsst/mops/removeSubsetsOMP.h"






namespace lsst {
namespace mops {





void populatePairsVectorFromFile(std::ifstream &pairsFile, 
                                 std::vector<std::set<unsigned int> > &pairsVector) {

    std::set<unsigned int>  tmpPair;
    std::string line;
    int tmpInt = -1;
    tmpPair.clear();
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
            tmpPair.insert(tmpInt);
        }
        if (tmpPair.size() < 2) {
            throw LSST_EXCEPT(InputFileFormatErrorException, "EE: CollapseTracklets: pairs in pairs file must be length >= 2!\n");
        }
        pairsVector.push_back(tmpPair);
        tmpPair.clear();
        
        line.clear();
        std::getline(pairsFile, line);   
    }
}




void writeTrackletsToOutFile(
    const std::vector<bool> & whichToKeep,
    const std::vector<std::set<unsigned int> > * tracklets, 
    std::ofstream &outFile){

    for (unsigned int i = 0; i < whichToKeep.size(); i++) {
        if (whichToKeep.at(i) == true) {
            std::set<unsigned int>::const_iterator indicesIter;
            for (indicesIter = tracklets->at(i).begin(); 
                 indicesIter != tracklets->at(i).end();
                 indicesIter++) {
                outFile << *indicesIter;
                outFile << " ";
            }
            outFile << std::endl;
        }
    }
}



  extern std::string curTime();

int removeSubsetsMain(int argc, char** argv) {

    /* read pairs from a file, do the removal, and write the output */
    std::string USAGE("USAGE: removeSubsets --inFile <input pairfile> --outFile <output pairfile>");
    std::ifstream inFile;
    std::ofstream outFile;
    std::vector <std::set<unsigned int> >* pairsVector = 
        new std::vector<std::set<unsigned int> >;
    std::vector <bool > outputVector;
    
    char* inFileName = NULL;
    char* outFileName = NULL;
    bool removeSubsets = true;
    bool shortCircuit = true;
    bool sortBeforeIntersect = false;
    
    static const struct option longOpts[] = {
        { "inFile", required_argument, NULL, 'i' },
        { "outFile", required_argument, NULL, 'o' },
        { "shortCircuit", required_argument, NULL, 'c'},
        { "sortBeforeIntersect", required_argument, NULL, 's'},
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
    };
    
    int longIndex = -1;
    const char* optString = "i:o:c:s:h";
    int opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
        case 'i':
            inFileName = optarg; 
            break;
            
        case 'o':
            outFileName = optarg ; 
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
    std::cout << "Short-circuit if possible:                   "<< boolToString(shortCircuit)
              << std::endl;
    std::cout << "Sort sets by size before set intersection:   "<< boolToString(sortBeforeIntersect)
              << std::endl;
    
    inFile.open(inFileName);
    outFile.open(outFileName);
    std::cout << "Reading infile, starting at " << curTime() << std::endl;
    populatePairsVectorFromFile(inFile, *pairsVector);
    std::cout << "Finished reading infile at " << curTime() << std::endl;
    
    /* do the actual work */
    
    
    SubsetRemover mySR;
    
    mySR.removeSubsetsPopulateOutputVector(pairsVector, 
                                           outputVector,
                                           shortCircuit, 
                                           sortBeforeIntersect);
    
    std::cout << "Finished filtering, writing output starting at " << curTime() << std::endl;
    writeTrackletsToOutFile(outputVector, pairsVector, outFile);
    outFile.close();

    delete pairsVector;
    if (outFile.good()) {
        std::cout << "Finished successfully finished at " << curTime() << std::endl;
        return 0;
    }
    else {
        std::cout << "ERROR writing/closing file. Finished at " << curTime() << std::endl;
        return 1;
    }
}



    
}} // close lsst::mops

int main(int argc, char** argv) {
    return lsst::mops::removeSubsetsMain(argc, argv);
}

// -*- LSST-C++ -*-

/*
 * jmyers 9/20/09
 * 
 * these were used for collapseTracklets back long ago, for testing
 * and to allow PanSTARRS to call collapseTracklets etc from a
 * shell. Though file I/O should no longer be done in C++ land in
 * production, I still use a lot of these tools from the command-line
 * for testing.
 */


#ifndef __MOPS_FILE_UTILS_H__
#define __MOPS_FILE_UTILS_H__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>


#include "MopsDetection.h"
#include "Exceptions.h"
#include "KDTree.h"
#include "Tracklet.h"


namespace lsst {
namespace mops {

/* if any index in pairs is >= detsSize, return false, else return true.
 * ALWAYS do this sanity check before anything else - particularly
 * setPhysicalParams vector or similar.
 */

bool isSane(unsigned int detsSize, const std::vector<Tracklet> *pairs);


void writeTrackletsToOutFile(const std::vector<Tracklet> * tracklets, std::ofstream &outFile);


void populateDetVectorFromFile(std::ifstream &detsFile, std::vector <MopsDetection> &myDets, const double &astromErr = 0.0);

void populatePairsVectorFromFile(std::ifstream &pairsFile,
				 std::vector <Tracklet> &pairsVector);

/* these overloaded versions are mainly for SWIG, since Python file objects != fstreams */
void writeTrackletsToOutFile(const std::vector<Tracklet> * tracklets, std::string outFileName);
    void populateDetVectorFromFile(std::string detsFileName, std::vector <MopsDetection> &myDets, const double &astromErr = 0.0);
void populatePairsVectorFromFile(std::string pairsFileName,
				 std::vector <Tracklet> &pairsVector);


}} // close namespace lsst::mops


#endif

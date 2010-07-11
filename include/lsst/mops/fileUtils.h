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


void populateDetVectorFromFile(std::ifstream &detsFile, std::vector <MopsDetection> &myDets);

void populatePairsVectorFromFile(std::ifstream &pairsFile,
				 std::vector <Tracklet> &pairsVector);

/* these overloaded versions are mainly for SWIG, since Python file objects != fstreams */
void writeTrackletsToOutFile(const std::vector<Tracklet> * tracklets, std::string outFileName);
void populateDetVectorFromFile(std::string detsFileName, std::vector <MopsDetection> &myDets);
void populatePairsVectorFromFile(std::string pairsFileName,
				 std::vector <Tracklet> &pairsVector);


}} // close namespace lsst::mops


#endif

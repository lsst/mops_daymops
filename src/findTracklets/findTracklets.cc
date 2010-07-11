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
 
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <getopt.h>
#include <set>
#include <map>
#include <string>
#include <sstream>
#include <math.h>

#include "lsst/mops/KDTree.h"
#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/daymops/findTracklets/findTracklets.h"

#define LEAF_NODE_SIZE 16

#define uint unsigned int

namespace lsst {
namespace mops {

/***
    Prototypes - these don't need to be seen by files which include findTracklets.h
***/

/*****************************************************************
 * Populate 2D vector of detections.
 * Each outer vector contains Detections of equal MJD.
 *****************************************************************/
void groupByImageTime(const std::vector<MopsDetection>&,
		      std::map<double, std::vector<MopsDetection> >&);


/******************************************************************
 * Take 2D vector of detections and create a map linking
 * each MJD vector to its MJD double value.
 ******************************************************************/
void generatePerImageTrees(const std::map<double, std::vector<MopsDetection> >&, 
                           std::map<double, KDTree<long int> >&);


/******************************************************************
 * Given a mapping of MJDs to KDTrees of PointAndValue pairs 
 * index by file line number index, generate tracklets for each 
 * query point within a distance determined by maxVelocity.
 ******************************************************************/
void getTracklets(std::vector<Tracklet>&, const std::map<double, KDTree<long int> >&, 
		  const std::vector<MopsDetection>&, double);







/*****************************************************************
 *The main function of this file.
 *****************************************************************/
std::vector<Tracklet> findTracklets(const std::vector<MopsDetection> &myDets, 
                                    double maxVelocity)
{
    //detection vectors, each vector of unique MJD
    std::map<double,  std::vector<MopsDetection> > detectionSets; 

    //vector of RA and dec pairs for later searching
    std::vector<double> queryPoints; 

    //link MJD to KDTree of unique MJD detections
    std::map<double, KDTree<long int> > myTreeMap; 

    groupByImageTime(myDets, 
                     detectionSets);
    
    generatePerImageTrees(detectionSets, myTreeMap);

    //get results
    std::vector<Tracklet> results;
    results.resize(0);
    getTracklets(results, myTreeMap, 
                 myDets, maxVelocity);

    return results;
}



/*****************************************************************
 * Populate 2D vector of detections.
 * Each outer vector contains Detections of equal MJD.
 *****************************************************************/
void groupByImageTime(const std::vector<MopsDetection> &myDets, 
		      std::map<double, std::vector<MopsDetection> > &detectionSets)
{
    // maps by default sort on their first parameter.
    // build a map from detection time to all dets at that time. we'll have a nice sorted
    // set of vectors.

    for (unsigned int i = 0; i < myDets.size(); i++) {
        detectionSets[myDets.at(i).getEpochMJD()].push_back(myDets.at(i));
    }
    
}


  

/******************************************************************
 * Take 2D vector of detections and create a map linking
 * each MJD vector to its MJD double value.
 ******************************************************************/
void generatePerImageTrees(const std::map<double, std::vector<MopsDetection> > &detectionSets, 
                           std::map<double, KDTree<long int> > &myTreeMap)
{

    // for each vector representing a single EpochMJD, created
    // a KDTree containing mapping detection RA, Decs -> detection IDs

    std::map<double, std::vector<MopsDetection> >::const_iterator imageIter;

    for(imageIter = detectionSets.begin(); imageIter != detectionSets.end(); imageIter++) {

        const std::vector<MopsDetection> *thisDetVec = &(imageIter->second);
        double thisEpoch = imageIter->first;
        std::vector<PointAndValue<long int> > vecPV;
        
        for(unsigned int j=0; j < thisDetVec->size(); j++) {

            PointAndValue<long int> tempPV;
            std::vector<double> pairRADec;
            
            pairRADec.push_back(convertToStandardDegrees(thisDetVec->at(j).getRA()));                    
            pairRADec.push_back(convertToStandardDegrees(thisDetVec->at(j).getDec()));
            
            tempPV.setPoint(pairRADec);
            tempPV.setValue(thisDetVec->at(j).getID());
            vecPV.push_back(tempPV);
        }
        
        KDTree<long int> tempKDTree;
        
        tempKDTree.buildFromData(vecPV, 2, LEAF_NODE_SIZE);
        myTreeMap.insert(std::make_pair(thisEpoch, tempKDTree) );
    }
}



/******************************************************************
 * Given a mapping of MJDs to KDTrees of PointAndValue pairs 
 * index by file line number index, generate tracklets for each 
 * query point within a distance determined by maxVelocity.
 ******************************************************************/
void getTracklets(std::vector<Tracklet> &resultsVec,  
		  const std::map<double, KDTree<long int> > &myTreeMap,
		  const std::vector<MopsDetection> &queryPoints,
		  double maxVelocity)
{
    // vectors of RADecRangeSearch parameters we search exclusively in RA, Dec;
    // the "otherDims" parameters sent to KDTree range search are empty.
    std::vector<double> otherDimsTolerances;
    std::vector<double> otherDimsPt;
    otherDimsTolerances.resize(0);
    otherDimsPt.resize(0);
    std::vector<GeometryType> myGeos;
    // we search RA, Dec only.
    myGeos.push_back(RA_DEGREES);
    myGeos.push_back(DEC_DEGREES);


    // iterate through list of collected Detections, as read from input
    // file, and search over them
    for(unsigned int i=0; i<queryPoints.size(); i++){
        
        const MopsDetection * curQuery = &(queryPoints.at(i));
        double queryRA = convertToStandardDegrees(curQuery->getRA());
        double queryDec = convertToStandardDegrees(curQuery->getDec());
        double queryMJD = curQuery->getEpochMJD();

        // iterate through each KDTree of detections, where each KDTree
        // represents a unique MJD
        std::map<double, KDTree<long int> >::const_iterator iter;
        for(iter = myTreeMap.begin(); iter != myTreeMap.end(); iter++) {

            double curMJD = iter->first;     //map key
            const KDTree<long int> *curTree = &(iter->second); //value associated with key

            //only consider this tree if it contains detections
            //that occurred after the current one
            if(curMJD > queryMJD){      
	  
                double maxDistance = (curMJD - queryMJD) * maxVelocity;
                std::vector<double> queryPt;
                queryPt.push_back(queryRA);
                queryPt.push_back(queryDec);
	
                std::vector<PointAndValue<long int> > queryResults;

                // do a rectangular search around this point. note that we 
                // use the haversine great-circle distance in the tree, so we are
                // sure that we get any object within maxDistance (and a few others)
                queryResults = curTree->RADecRangeSearch(queryPt, maxDistance,
                                                         otherDimsPt, otherDimsTolerances,
                                                         myGeos);
                
                // filter the results, getting the items which are actually within
                // the circle we are searching, not the rectangle enclosing it.
                std::vector<long int> closeEnoughResults;
                    
                for (unsigned int ii = 0; ii < queryResults.size(); ii++) {
                    PointAndValue<long int> * curResult = &(queryResults.at(ii));
                    double properDistance =  angularDistanceRADec_deg(queryRA, 
                                                                      queryDec, 
                                                                      curResult->getPoint().at(0),
                                                                      curResult->getPoint().at(1));
                    if (properDistance <= maxDistance) {
                        closeEnoughResults.push_back(curResult->getValue());
                    }
                }
                
                // collect results for each query point's results for each MJD
                Tracklet newTracklet;
                newTracklet.indices.insert(curQuery->getID());               
                
                for(unsigned int j=0; j < closeEnoughResults.size(); j++){
                    newTracklet.indices.insert(closeEnoughResults.at(j));
                }

                if (newTracklet.indices.size() >= 2) {
                    resultsVec.push_back(newTracklet);
                }
                queryResults.clear();
                queryPt.clear();
            }
        }
    }
}


}} // close lsst::mops

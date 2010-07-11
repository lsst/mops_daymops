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
 
/* File: detectionProximity.cc
 * Author: Matthew Cleveland
 * Purpose: 
 */

#include "lsst/mops/daymops/detectionProximity/detectionProximity.h"


namespace lsst {
    namespace mops {

// prototypes not to be seen outside this file

void buildKDTree(const std::vector<MopsDetection>, KDTree<unsigned int>&);

std::vector<std::pair <unsigned int, unsigned int> > getProximity(const std::vector<MopsDetection>& queryPoints,
								  const KDTree<unsigned int>& searchTree,
								  double maxDist,
								  double maxTime);


std::vector<std::pair <unsigned int, unsigned int> > detectionProximity(
    const std::vector<MopsDetection>& queryPoints,
    const std::vector<MopsDetection>& dataPoints,
    double distanceThreshold,
    double timeThreshold)
{
    std::vector<std::pair <unsigned int, unsigned int> > results;
    
    if(queryPoints.size() > 0 && dataPoints.size() > 0){
        
        //build KDTrees from detection vectors
        KDTree<unsigned int> dataTree;
        buildKDTree(dataPoints, dataTree);
        
        //get results
        results = getProximity(queryPoints, dataTree, distanceThreshold,
                               timeThreshold);
    }

    return results;
}



/**********************************************************************
 * Populate a KDTree 'tree' from the Detections in vector 'points'
 ***********************************************************************/
void buildKDTree(const std::vector<MopsDetection> points, KDTree<unsigned int> &tree)
{

  if(points.size() > 0){

    std::vector<PointAndValue<unsigned int> > vecPV;

    for(unsigned int i=0; i < points.size(); i++){
    
	PointAndValue<unsigned int> tempPV;
	std::vector<double> pairRADec;
	
	pairRADec.push_back(convertToStandardDegrees(points.at(i).getRA()));                    
	pairRADec.push_back(convertToStandardDegrees(points.at(i).getDec()));
        pairRADec.push_back(points.at(i).getEpochMJD());
	
	tempPV.setPoint(pairRADec);
	tempPV.setValue(i);
	vecPV.push_back(tempPV);
    }
    
    tree.buildFromData(vecPV, 4, 100);
  }
}
    
    

/*
 *
 */
std::vector<std::pair <unsigned int, unsigned int> > getProximity(const std::vector<MopsDetection>& queryPoints,
								  const KDTree<unsigned int>& searchTree,
								  double maxDist,
								  double maxTime)
{
  std::vector<std::pair <unsigned int, unsigned int> > pairs;

  std::vector<GeometryType> myGeos;
  myGeos.push_back(RA_DEGREES); //RA
  myGeos.push_back(DEC_DEGREES); //Dec
  myGeos.push_back(EUCLIDEAN); //time

  
  for(unsigned int i=0; i<queryPoints.size(); i++){

      std::vector<double> RADecQueryPt;
      std::vector<double> otherDimsPt;
      std::vector<double> otherDimsTolerances;
      std::vector<PointAndValue<unsigned int> > queryResults;

      RADecQueryPt.push_back(convertToStandardDegrees(queryPoints.at(i).getRA()));
      RADecQueryPt.push_back(convertToStandardDegrees(queryPoints.at(i).getDec()));
    
      otherDimsPt.push_back(queryPoints.at(i).getEpochMJD());
      
      otherDimsTolerances.push_back(maxTime);
      
      queryResults = searchTree.RADecRangeSearch(RADecQueryPt, maxDist, 
                                                 otherDimsPt, otherDimsTolerances, 
                                                 myGeos);
      
      for(unsigned int j=0; j<queryResults.size(); j++){
          pairs.push_back(std::make_pair(i, queryResults.at(j).getValue()));
      }
  }
  
  
  return pairs;
}


    }} // close lsst::mops

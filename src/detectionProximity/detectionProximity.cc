// -*- LSST-C++ -*-
/* File: detectionProximity.cc
 * Author: Matthew Cleveland
 * Purpose: 
 */

#include "detectionProximity.h"



// prototypes not to be seen outside this file

void buildKDTree(const std::vector<Detection>, KDTree::KDTree<unsigned int>&);

std::vector<std::pair<unsigned int, unsigned int> > 
getProximity(const std::vector<Detection>& queryPoints,
	     const KDTree::KDTree<unsigned int>& searchTree, double,
	     double, double);



std::vector<std::pair <unsigned int, unsigned int> > detectionProximity(const std::vector<Detection>& queryPoints,
								   const std::vector<Detection>& dataPoints,
								   double distanceThreshold,
								   double brightnessThreshold,
								   double timeThreshold)
{
    std::vector<std::pair <unsigned int, unsigned int> > results;
    
    if(queryPoints.size() > 0 && dataPoints.size() > 0){
        
        //build KDTrees from detection vectors
        KDTree::KDTree<unsigned int> dataTree;
        buildKDTree(dataPoints, dataTree);
        
        //get results
        results = getProximity(queryPoints, dataTree, distanceThreshold,
                               brightnessThreshold, timeThreshold);
    }

    return results;
}



/**********************************************************************
 * Populate a KDTree 'tree' from the Detections in vector 'points'
 ***********************************************************************/
void buildKDTree(const std::vector<Detection> points, KDTree::KDTree<unsigned int> &tree)
{

  if(points.size() > 0){

    std::vector<KDTree::PointAndValue<unsigned int> > vecPV;

    for(unsigned int i=0; i < points.size(); i++){
    
	KDTree::PointAndValue<unsigned int> tempPV;
	std::vector<double> pairRADec;
	
	pairRADec.push_back(KDTree::Common::convertToStandardDegrees(points.at(i).getRA()));                    
	pairRADec.push_back(KDTree::Common::convertToStandardDegrees(points.at(i).getDec()));
        pairRADec.push_back(points.at(i).getMag());
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
std::vector<std::pair <unsigned int, unsigned int> > getProximity(const std::vector<Detection>& queryPoints,
								  const KDTree::KDTree<unsigned int>& searchTree,
								  double maxDist, double maxBrightness,
								  double maxTime)
{
  std::vector<std::pair <unsigned int, unsigned int> > pairs;

  std::vector<KDTree::Common::GeometryType> myGeos;
  myGeos.push_back(KDTree::Common::RA_DEGREES); //RA
  myGeos.push_back(KDTree::Common::DEC_DEGREES); //Dec
  myGeos.push_back(KDTree::Common::EUCLIDEAN); //magnitude
  myGeos.push_back(KDTree::Common::EUCLIDEAN); //time

  
  for(unsigned int i=0; i<queryPoints.size(); i++){

      std::vector<double> RADecQueryPt;
      std::vector<double> otherDimsPt;
      std::vector<double> otherDimsTolerances;
      std::vector<KDTree::PointAndValue<unsigned int> > queryResults;

      RADecQueryPt.push_back(KDTree::Common::convertToStandardDegrees(queryPoints.at(i).getRA()));
      RADecQueryPt.push_back(KDTree::Common::convertToStandardDegrees(queryPoints.at(i).getDec()));
    
      otherDimsPt.push_back(queryPoints.at(i).getMag());
      otherDimsPt.push_back(queryPoints.at(i).getEpochMJD());
      
      otherDimsTolerances.push_back(maxBrightness);
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

// -*- LSST-C++ -*-
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


namespace lsst {
    namespace mops {

/***
    Prototypes - these don't need to be seen by files which include findTracklets.h
***/

/*****************************************************************
 * Populate 2D vector of detections.
 * Each outer vector contains Detections of equal MJD.
 *****************************************************************/
void generateMJDTrees(const std::vector<MopsDetection>*,
		      std::vector< std::vector<MopsDetection> >*);


/******************************************************************
 * Take 2D vector of detections and create a map linking
 * each MJD vector to its MJD double value.
 ******************************************************************/
void generateTreeMap(std::vector< std::vector<MopsDetection> >*, 
		     std::map<double, KDTree<int> >*);


/******************************************************************
 * Given a mapping of MJDs to KDTrees of PointAndValue pairs 
 * index by file line number index, generate tracklets for each 
 * query point within a distance determined by maxVelocity.
 ******************************************************************/
void getTracklets(std::vector<int>*, std::map<double, KDTree<int> >*, 
		  const std::vector<MopsDetection>*, double);




/**************************************************************************
 * Find the tracklets within minVelocity distance of the query points
 * and then remove them from the set of tracklets within maxVelocity
 * distance of query points.
 **************************************************************************/
void prune(std::vector<int>*, std::map<double, KDTree<int> >*,
	   const std::vector<MopsDetection>*, double);


/****************************************************************
 * Find the position of a double in a vector of doubles.  Return
 * position index or -1 if not found.
 ****************************************************************/
int vectorPosition(const std::vector<double>*, double);


/**************************************************************
 * Determine all unique epochs found in a set of detections.
 * Populate "epochs" vector with these unique values.
 **************************************************************/
void getEpochs(const std::vector<MopsDetection>*, std::vector<double>*);


/*****************************************************************
 *The main function of this file.
 *****************************************************************/
std::vector<Tracklet> findTracklets(const std::vector<MopsDetection> &myDets, 
                                    double maxVelocity, double minVelocity)
{
    //detection vectors, each vector of unique MJD
    std::vector< std::vector<MopsDetection> > detectionSets; 

    //vector of RA and dec pairs for later searching
    std::vector<double> queryPoints; 

    //link MJD to KDTree of unique MJD detections
    std::map<double, KDTree<int> > myTreeMap; 

    generateMJDTrees(&myDets, &detectionSets);

    generateTreeMap(&detectionSets, &myTreeMap);

    //get results
    std::vector<int> results;
    results.resize(0);
    getTracklets(&results, &myTreeMap, 
                 &myDets, maxVelocity);

    //prune results if minimum velocity is specified
    if(minVelocity > 0){
        prune(&results, &myTreeMap, &myDets, minVelocity);
    }

    //Create Tracklet vector from results vector and return it.
    std::vector<Tracklet> tracklets;
  
    for(unsigned int i=0; i<results.size(); i++){

        Tracklet tempTracklet;
        tempTracklet.indices.insert(results.at(i));
        i++;
        tempTracklet.indices.insert(results.at(i));
        tracklets.push_back(tempTracklet);
    }
  
    return tracklets;
}



/*****************************************************************
 * Populate 2D vector of detections.
 * Each outer vector contains Detections of equal MJD.
 *****************************************************************/
void generateMJDTrees(const std::vector<MopsDetection> *myDets, 
		      std::vector< std::vector<MopsDetection> > *detectionSets)
{
    if(myDets->size() > 0){

        std::vector<double> epochs;
        epochs.resize(0); //initialize
        getEpochs(myDets, &epochs);
    
        double epochTime;

        //initialize vector to number of different epochs found
        detectionSets->resize(epochs.size());
    
        int numDetections = myDets->size();
        int insertIndex;

        MopsDetection tempD;
    
        //populate 2D vector, ordered by EpochMJD, from 1D 
        //vector myDets
        for(int i=0; i < numDetections; i++){
            epochTime = myDets->at(i).getEpochMJD();

            //get the current detection and set its file index value
            tempD = myDets->at(i);
            tempD.setFileIndex(i);

            insertIndex = vectorPosition(&epochs, epochTime);

            detectionSets->at(insertIndex).push_back(tempD);
        }
    }
}


  

/******************************************************************
 * Take 2D vector of detections and create a map linking
 * each MJD vector to its MJD double value.
 ******************************************************************/
void generateTreeMap(std::vector< std::vector<MopsDetection> > *detectionSets, 
		     std::map<double, KDTree<int> > *myTreeMap)
{
    // for each vector representing a single EpochMJD, created
    // a KDTree containing its line number in the input file 
    // and PointAndValue pair
    for(unsigned int i=0; i < detectionSets->size(); i++){
        int thisTreeSize = detectionSets->at(i).size();
    
        if(thisTreeSize > 0){
            std::vector<MopsDetection> thisDetVec = detectionSets->at(i);

            if(thisDetVec.size() > 0){
                double thisEpoch = thisDetVec.at(0).getEpochMJD();
      
                std::vector<PointAndValue<int> > vecPV;
	
                for(int j=0; j < thisTreeSize; j++){
                    PointAndValue<int> tempPV;
                    std::vector<double> pairRADec;

                    pairRADec.push_back(convertToStandardDegrees(thisDetVec.at(j).getRA()));                    
                    pairRADec.push_back(convertToStandardDegrees(thisDetVec.at(j).getDec()));
	  
                    tempPV.setPoint(pairRADec);
                    tempPV.setValue(thisDetVec.at(j).getFileIndex());
                    vecPV.push_back(tempPV);
                }
                KDTree<int> tempKDTree;

                tempKDTree.buildFromData(vecPV, 2, 100);
                myTreeMap->insert(std::make_pair(thisEpoch, tempKDTree) );
            }
        }
    }
}



/******************************************************************
 * Given a mapping of MJDs to KDTrees of PointAndValue pairs 
 * index by file line number index, generate tracklets for each 
 * query point within a distance determined by maxVelocity.
 ******************************************************************/
void getTracklets(std::vector<int> *resultsVec,  
		  std::map<double, KDTree<int> > *myTreeMap,
		  const std::vector<MopsDetection> *queryPoints,
		  double maxVelocity)
{
    // vectors of RADecRangeSearch parameters
    // we search exclusively in RA, Dec, so the 'other dimensions' are all empty
    std::vector<double> otherDimsTolerances;
    std::vector<double> otherDimsPt;
    otherDimsTolerances.resize(0);
    otherDimsPt.resize(0);
    double maxDistance;
    std::vector<GeometryType> myGeos;
    // we search RA, Dec only.
    myGeos.push_back(RA_DEGREES);
    myGeos.push_back(DEC_DEGREES);

    // hyperRectangleSearch result container
    std::vector<PointAndValue<int> > queryResults;
    std::map<double, KDTree<int> >::iterator iter;

    // loop variables
    int leftIndex, rightIndex, count=0;
    double treeMJD, tempDec, tempRA, tempMJD;
    KDTree<int> tempKDTree;
    std::vector<double> queryPt;
    PointAndValue<int> tempPV;
    MopsDetection tempD;

    // iterate through list of collected Detections, as read from input
    // file, and search over them
    for(unsigned int i=0; i<queryPoints->size(); i++){

        // iterate through each KDTree of detections, where each KDTree
        // represents a unique MJD
        for(iter = myTreeMap->begin(); iter != myTreeMap->end(); ++iter) {

            tempD = queryPoints->at(i);
            tempRA = convertToStandardDegrees(tempD.getRA());
            tempDec = convertToStandardDegrees(tempD.getDec());
            tempMJD = tempD.getEpochMJD();
	
            treeMJD = iter->first;     //map key
            tempKDTree = iter->second; //value associated with key

            //only consider this tree if it contains detections
            //that occurred after the current one
            if(treeMJD > tempMJD){      
	  
                maxDistance = (treeMJD - tempMJD) * maxVelocity;

                queryPt.push_back(tempRA);
                queryPt.push_back(tempDec);

                leftIndex = i; //index of query point in input file
	
                queryResults = tempKDTree.RADecRangeSearch(queryPt, maxDistance,
                                                           otherDimsPt, otherDimsTolerances,
                                                           myGeos);

                // collect results for each query point's results for each MJD
                for(unsigned int j=0; j<queryResults.size(); j++){
	  
                    rightIndex = queryResults.at(j).getValue();                    
                    resultsVec->push_back(leftIndex);
                    resultsVec->push_back(rightIndex);
                }
	
                count = i; /* count is used so that when iterating through query points, those
                            * that are of a MJD less than the current one are not considered */
            }
            queryPt.clear();
            queryResults.clear();
        }
    }
}




/**************************************************************************
 * Find the tracklets within minVelocity distance of the query points
 * and then remove them from the set of tracklets within maxVelocity
 * distance of query points.
 **************************************************************************/
void prune(std::vector<int> *resultsVec, 
	   std::map<double, KDTree<int> > *myTreeMap,
	   const std::vector<MopsDetection> *queryPoints, double minVelocity)
{
    std::set <std::vector<int> > maxSet, minSet, diffSet;
  
    std::vector<int> insertTemp;

    for(unsigned int i=0; i < resultsVec->size();i++){
        insertTemp.push_back(resultsVec->at(i));
        i++;
        insertTemp.push_back(resultsVec->at(i));
        maxSet.insert(insertTemp);
        insertTemp.clear();
    }

    std::vector<int> minResults;
    minResults.resize(0);
    getTracklets(&minResults, myTreeMap, queryPoints, minVelocity);
  
    for(unsigned int i=0; i < minResults.size(); i++){
        insertTemp.push_back(minResults.at(i));
        i++;
        insertTemp.push_back(minResults.at(i));
        minSet.insert(insertTemp);
        insertTemp.clear();
    }

    std::insert_iterator <std::set <std::vector<int> > > ii(diffSet, diffSet.begin());
  
    std::set_difference(maxSet.begin(), maxSet.end(), minSet.begin(), minSet.end(), ii);

    resultsVec->clear();
  
    std::set <std::vector<int> >::iterator iter;
    for(iter = diffSet.begin(); iter != diffSet.end(); iter++){
        resultsVec->push_back(iter->at(0));
        resultsVec->push_back(iter->at(1));
    }
}





/**************************************************************
 * Determine all unique epochs found in a set of detections.
 * Populate "epochs" vector with these unique values.
 **************************************************************/
void getEpochs(const std::vector<MopsDetection> *myDets, 
	       std::vector<double> *epochs)
{
    double tempEpoch;
    int vPos;

    for(unsigned int i=0; i<myDets->size(); i++){
        tempEpoch=myDets->at(i).getEpochMJD();
    
        vPos = vectorPosition(epochs, tempEpoch);

        if(vPos == -1){
            epochs->push_back(tempEpoch);
        }
    }
}




/****************************************************************
 * Find the position of a double in a vector of doubles.  Return
 * position index or -1 if not found.
 ****************************************************************/
int vectorPosition(const std::vector<double>* lookUp, double val)
{
    std::vector<double>::iterator result;
  
    unsigned int position = std::find(lookUp->begin(), lookUp->end(), val) - lookUp->begin();

    if(position < lookUp->size()){
        return position;
    }
    else{
        return -1;
    }
}





    }} // close lsst::mops

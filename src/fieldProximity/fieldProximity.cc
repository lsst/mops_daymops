// -*- LSST-C++ -*-
#include <stdlib.h>
#include <stdio.h>
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
#include "lsst/mops/daymops/fieldProximity/fieldProximity.h"


namespace lsst {
    namespace mops {

// private declarations

std::vector<FieldProximityTrack> buildEphemeridesVec(const std::vector<FieldProximityTrack>,
                                     std::vector<std::string>*,
                                     const std::vector<double>);

void generateMJDVecs(const std::vector<FieldProximityTrack>, 
                     std::vector<double>,
                     std::vector<std::vector<FieldProximityTrack> >*);

void generateTreeMap(std::vector< std::vector<FieldProximityTrack> >*, 
		     std::map<double, KDTree<int> >*);

void getProximity(std::vector<int>*,
                  std::map<double, KDTree<int> >*,
                  const std::vector<Field>);

int vecPosition(const std::vector<std::string>*, std::string);

int vecPosition(const std::vector<double>*, double);

std::vector<double> getFieldTimes(std::vector<Field>, std::vector<std::string>*);

std::vector<std::pair <unsigned int, unsigned int> > makePair(std::vector<int>);


/**************************************************
 *
 *                 main function
 *
 **************************************************/
std::vector<std::pair <unsigned int, unsigned int> > 
fieldProximity(std::vector<FieldProximityTrack> allTracks,
               std::vector<Field> queryFields,
               double distThresh) 
{
    //return value
    std::vector<std::pair <unsigned int, unsigned int> > pairsResult;

    if(queryFields.size() > 0 && allTracks.size() > 0){
        std::vector<std::string> trackIDs, fieldIDs;
        
        //a vector of FPT vectors, each outer vector contains
        //a vector of FPT of a specific time
        std::vector<std::vector<FieldProximityTrack> > tracksByTime;
        
        //link MJD to KDTree of unique MJD detections
        std::map<double, KDTree<int> > myTreeMap; 
        
        std::vector<double> fieldTimes = getFieldTimes(queryFields, &fieldIDs);
        
        std::vector<FieldProximityTrack> ephemeridesTracks = buildEphemeridesVec(allTracks,
                                                                                 &trackIDs,
                                                                                 fieldTimes);
        
        generateMJDVecs(ephemeridesTracks, fieldTimes, &tracksByTime);
        
        generateTreeMap(&tracksByTime, &myTreeMap);
        
        //get results
        std::vector<int> results;
      
        getProximity(&results, &myTreeMap, queryFields);
        
        //format results into desired, pair, output
        pairsResult = makePair(results);
    }
    
    return pairsResult;
}





/***************************************************
 *
 ***************************************************/
std::vector<double> getFieldTimes(std::vector<Field> fields,
                                  std::vector<std::string>* fieldIDs)
{
    std::vector<double> results;

    for(unsigned int i=0; i < fields.size(); i++){
        double time = fields.at(i).getEpochMJD();

        //keep track of all field times to be searched
        int vecPos = vecPosition(&results, time);
        
        if(vecPos == -1){
            results.push_back(time);
        }

        fieldIDs->push_back(fields.at(i).getFieldID());
                          
    }

    return results;
}



/***************************************************
 *
 ***************************************************/
std::vector<FieldProximityTrack> buildEphemeridesVec(const std::vector<FieldProximityTrack> tracks,
                                       std::vector<std::string> *trackIDs,
                                       const std::vector<double> fieldTimes)
{
    std::vector<FieldProximityTrack> finalTracks;
    
    double pt1RA, pt1Dec, pt1Time;
    double pt2RA, pt2Dec, pt2Time;
    
    std::string trackID;
    
    //find track link points and calculate for each track
    for(unsigned int i=0; i < tracks.size(); i++){

        trackID = tracks.at(i).getID();

        //keep a list of input tracks
        trackIDs->push_back(trackID);

        std::vector<FieldProximityPoint> points = tracks.at(i).getPoints();

        //first point
        pt1RA = points.at(0).getRA();
        pt1Dec = points.at(0).getDec();
        pt1Time = points.at(0).getEpochMJD();

        pt2RA = points.at(1).getRA();
        pt2Dec = points.at(1).getDec();
        pt2Time = points.at(1).getEpochMJD();
    
        //generate the ephemerides points along the line described
        //by pts 1 and 2 for each possible time
        for(unsigned int j=0; j < fieldTimes.size(); j++){
            
            double ephRA, ephDec;

            //for projecting forward in time
            if(fieldTimes.at(j) > pt2Time){
                ephRA = pt2RA + (pt2RA - pt1RA) * (fieldTimes.at(j) - pt2Time);
                ephDec = pt2Dec + (pt2Dec - pt1Dec) * (fieldTimes.at(j) - pt2Time);
            }
            //for projecting backwards in time
            else if(fieldTimes.at(j) < pt1Time){
                ephRA = pt1RA - (pt2RA - pt1RA) * (pt1Time - fieldTimes.at(j));
                ephDec = pt1RA - (pt2Dec - pt1Dec) * (pt1Time - fieldTimes.at(j));
            }
            //for ephemera in between pt1 and pt2 times
            else{
                double timeFactor = (pt2Time - fieldTimes.at(j))/(pt2Time - pt1Time);
                ephRA = pt1RA + ((pt2RA - pt1RA) * timeFactor);
                ephDec = pt1Dec + ((pt2Dec - pt1Dec) * timeFactor);
            }

            FieldProximityTrack f;
            std::vector<FieldProximityPoint> addPoints;
            FieldProximityPoint pt;
            
            pt.setRA(convertToStandardDegrees(ephRA));
            pt.setDec(convertToStandardDegrees(ephDec));
            pt.setEpochMJD(fieldTimes.at(j)); //i ???
            
            addPoints.push_back(pt);
            
            f.setPoints(addPoints);
            
            std::string s;
            std::stringstream out;
            out << i;
            f.setID(out.str());
            
            finalTracks.push_back(f);
        }
    }
    
    return finalTracks;
}



/*****************************************************************
 * Create a vector of Track vectors.  Each inner vector will 
 * contain all tracks of a specific MJD.
 *****************************************************************/
void generateMJDVecs(const std::vector<FieldProximityTrack> allTracks, 
                     std::vector<double> fieldTimes,
                     std::vector<std::vector<FieldProximityTrack> > *tracksByTime)
{
    double thisTime;
    int vecPos;

    tracksByTime->resize(fieldTimes.size());

    for(unsigned int i=0; i < allTracks.size(); i++){
        
        thisTime = allTracks.at(i).getPoints().at(0).getEpochMJD();
        
        vecPos = vecPosition(&fieldTimes, thisTime);

        if(vecPos != -1){
            FieldProximityTrack f = allTracks.at(i);
            tracksByTime->at(vecPos).push_back(f);            
        }
    }
}


  

/******************************************************************
 * Take 2D vector of detections and create a map linking
 * each MJD vector to its MJD double value.
 ******************************************************************/
void generateTreeMap(std::vector< std::vector<FieldProximityTrack> > *trackSets, 
		     std::map<double, KDTree<int> > *myTreeMap)
{
    // for each vector representing a single EpochMJD, created
    // a KDTree containing its line number in the input file 
    // and PointAndValue pair
    for(unsigned int i=0; i < trackSets->size(); i++){

        int thisTreeSize = trackSets->at(i).size();

        if(thisTreeSize > 0){
            
            std::vector<FieldProximityTrack> thisTrackVec = trackSets->at(i);

            if(thisTrackVec.size() > 0){

                FieldProximityTrack t;
                t = thisTrackVec.at(0);
                std::vector<FieldProximityPoint> pts;
                pts = t.getPoints();

                double thisEpoch = pts.at(0).getEpochMJD();
      
                std::vector<PointAndValue<int> > vecPV;
	
                for(int j=0; j < thisTreeSize; j++){

                    PointAndValue<int> tempPV;
                    std::vector<double> pairRADec;

                    t = thisTrackVec.at(j);
                    pts = t.getPoints();

                    pairRADec.push_back(convertToStandardDegrees(pts.at(0).getRA()));
                    pairRADec.push_back(convertToStandardDegrees(pts.at(0).getDec()));
                    
                    tempPV.setPoint(pairRADec);
                    tempPV.setValue(atoi(t.getID().c_str()));
                    vecPV.push_back(tempPV);

                }
                KDTree<int> tempKDTree(vecPV, 2, 100);
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
void getProximity(std::vector<int> *resultsVec,  
		  std::map<double, KDTree<int> > *myTreeMap,
		  const std::vector<Field> queryPoints)
{
    // vectors of hyperRectangleSearch parameters
    std::vector<double> tolerances;
  
    std::vector<GeometryType> myGeos;
    myGeos.push_back(CIRCULAR_DEGREES);
    myGeos.push_back(CIRCULAR_DEGREES);

    // hyperRectangleSearch result container
    std::vector<PointAndValue<int> > queryResults;
    std::map<double, KDTree<int> >::iterator iter;

    //temps for querying
    Field tempField;
    double tempRA, tempDec, tempMJD, tempRadius, treeMJD;
    int tempFileIndex, leftIndex, rightIndex;
    KDTree<int> tempKDTree;
    std::vector<double> queryPt;

    //iterate through all fields collected, querying them against
    //the appropriate KDTrees
    for(unsigned int i=0; i<queryPoints.size(); i++){

        // iterate through each KDTree of tracks, where each KDTree
        // represents a unique MJD
        for(iter = myTreeMap->begin(); iter != myTreeMap->end(); ++iter) {

            tempField = queryPoints.at(i);

            tempRA = convertToStandardDegrees(tempField.getRA());
            tempDec = convertToStandardDegrees(tempField.getDec());
            tempMJD = tempField.getEpochMJD();
            tempFileIndex = tempField.getFileIndex();
            tempRadius = tempField.getRadius();

            treeMJD = iter->first;     //map key
            tempKDTree = iter->second; //value associated with key

            //only consider this tree if it contains fields
            //that occurred at the considered time
            if(treeMJD == tempMJD){
	  
                tolerances.push_back(tempRadius);
                tolerances.push_back(tempRadius);

                queryPt.push_back(tempRA);
                queryPt.push_back(tempDec);

                leftIndex = i; //index of query Field considered
	
                queryResults = tempKDTree.hyperRectangleSearch(queryPt, tolerances, myGeos);

                // collect results for each query point's results for each MJD
                for(unsigned int j=0; j<queryResults.size(); j++){
	  
                    rightIndex = queryResults.at(j).getValue();//
                    
                    std::vector<double> tempv = queryResults.at(j).getPoint();
                    
                    resultsVec->push_back(leftIndex);
                    resultsVec->push_back(rightIndex);
                }
        
                tolerances.clear();
                queryPt.clear();
                queryResults.clear();
            }
        }
    }
}




/****************************************************************
 * Find the position of a double in a vector of strings.  Return
 * position index or -1 if not found.
 ****************************************************************/
int vecPosition(const std::vector<std::string>* lookUp, std::string val)
{
    std::vector<std::string>::iterator result;
  
    unsigned int position = std::find(lookUp->begin(), lookUp->end(), val) - lookUp->begin();

    if(position < lookUp->size()){
        return position;
    }
    else{
        return -1;
    }
}

/***************************************************
 * Find the position of a double in a vector of doubles.  Return
 * position index or -1 if not found.
 ***************************************************/
int vecPosition(const std::vector<double>* lookUp, double val)
{
    std::vector<std::string>::iterator result;
  
    unsigned int position = std::find(lookUp->begin(), lookUp->end(), val) - lookUp->begin();

    if(position < lookUp->size()){
        return position;
    }
    else{
        return -1;
    }
}




/***************************************************
 *
 ***************************************************/
std::vector<std::pair <unsigned int, unsigned int> > makePair(std::vector<int> results)
{
    std::vector<std::pair <unsigned int, unsigned int> > returnVal;

    for(unsigned int i=0; i<results.size(); i++){
        std::pair <unsigned int, unsigned int> temp;
        temp.first = results.at(i);
        i++;
        temp.second = results.at(i);
        returnVal.push_back(temp);
    }

    return returnVal;
}


}} // close lsst::mops 

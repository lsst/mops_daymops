// -*- LSST-C++ -*-
/* File: orbitProximity.cc
 * Author: Matthew Cleveland
 * Purpose: 
 */

#include "../KDTree.h"
#include "orbitProximity.h"


#define DEBUG false

// "internal" declarations

KDTree::KDTree<unsigned int> buildKDTree(const std::vector<Orbit>);



/* use const by reference to avoid copying lots of memory */
std::vector<std::pair <unsigned int, unsigned int>  > 
getResults(const KDTree::KDTree<unsigned int>&,
           const std::vector<Orbit>&,
           double, double, double,
           double, double, double);










std::vector<std::pair <unsigned int, unsigned int> > 
orbitProximity(std::vector<Orbit> dataOrbits, 
               std::vector<Orbit> queryOrbits,
               double perihelionTolerance,
               double eccentricityTolerance,
               double inclinationTolerance,
               double perihelionArgTolerance,
               double longitudeArgTolerance,
               double perihelionTimeTolerance)
{
    //build KDTrees from Orbit vectors
    KDTree::KDTree<unsigned int> dataTree;
    dataTree = buildKDTree(dataOrbits);
    
    if (DEBUG) {
        std::cerr << "dataOrbits size: " << dataOrbits.size() << std::endl;
        std::cerr << "queryOrbits size: " << queryOrbits.size() << std::endl;
    }


    //get results
    std::vector<std::pair <unsigned int, unsigned int > >
        results = getResults(dataTree, queryOrbits,
                             perihelionTolerance, 
                             eccentricityTolerance, 
                             inclinationTolerance, 
                             perihelionArgTolerance, 
                             longitudeArgTolerance, 
                             perihelionTimeTolerance);
    
    return results;
}







/**********************************************************************
 * Populate a KDTree 'tree' from the Orbits in vector 'orbits'
 ***********************************************************************/
KDTree::KDTree<unsigned int> buildKDTree(const std::vector<Orbit> orbits)
{
    KDTree::KDTree<unsigned int> toReturn;
    
    std::vector<KDTree::PointAndValue<unsigned int> > vecPV;
    
    const int orbitDimensions = 6;
    
    for(unsigned int i=0; i < orbits.size(); i++){
        
	KDTree::PointAndValue<unsigned int> tempPV;
	std::vector<double> orbitDims;

	orbitDims.push_back(orbits.at(i).getPerihelion());
	orbitDims.push_back(orbits.at(i).getEccentricity());
	orbitDims.push_back(orbits.at(i).getInclination());
	orbitDims.push_back(orbits.at(i).getPerihelionArg());
	orbitDims.push_back(orbits.at(i).getLongitude());
	orbitDims.push_back(orbits.at(i).getPerihelionTime());
	
	tempPV.setPoint(orbitDims);
	tempPV.setValue(i);
	vecPV.push_back(tempPV);
    }
    
    toReturn.buildFromData(vecPV, orbitDimensions, 100);
    return toReturn;
}



/**************************************************
 *
 **************************************************/
std::vector<std::pair<unsigned int, unsigned int> > 
getResults(const KDTree::KDTree<unsigned int> &searchTree,
           const std::vector<Orbit> &queryPoints,
           double maxPerihelion,
           double maxEccentricity,
           double maxInclination,
           double maxPerihelionArg,
           double maxLongitude,
           double maxPerihelionTime)
{
    std::vector<std::pair<unsigned int, unsigned int> > results;
    std::vector<KDTree::PointAndValue<unsigned int> > queryResults;
    std::vector<double> queryPt;

    std::vector<KDTree::Common::GeometryType> myGeos;

    myGeos.push_back(KDTree::Common::CIRCULAR_DEGREES); //perihelion
    myGeos.push_back(KDTree::Common::EUCLIDEAN);        //eccentricity 
    myGeos.push_back(KDTree::Common::CIRCULAR_DEGREES); //inclination
    myGeos.push_back(KDTree::Common::CIRCULAR_DEGREES); //arg. of perihelion
    myGeos.push_back(KDTree::Common::CIRCULAR_DEGREES); //longitude
    myGeos.push_back(KDTree::Common::EUCLIDEAN);        //time of perihelion

    std::vector<double> tolerances;
    tolerances.push_back(maxPerihelion);
    tolerances.push_back(maxEccentricity);
    tolerances.push_back(maxInclination);
    tolerances.push_back(maxPerihelionArg);
    tolerances.push_back(maxLongitude);
    tolerances.push_back(maxPerihelionTime);

    for(unsigned int i=0; i<queryPoints.size(); i++){

        // degrees - perihelion (periapsis)
        queryPt.push_back(
            KDTree::Common::convertToStandardDegrees(queryPoints.at(i).getPerihelion()));
        //euclidean - eccentricity
        queryPt.push_back(queryPoints.at(i).getEccentricity());
        // degrees - inclination
        queryPt.push_back(
            KDTree::Common::convertToStandardDegrees(queryPoints.at(i).getInclination()));
        //degrees - arg. of perihelion (periapsis)
        queryPt.push_back(
            KDTree::Common::convertToStandardDegrees(queryPoints.at(i).getPerihelionArg()));
        //degrees - longitude
        queryPt.push_back(
            KDTree::Common::convertToStandardDegrees(queryPoints.at(i).getLongitude()));
        //euclidean - perihelion time (time of periapsis)
        queryPt.push_back(queryPoints.at(i).getPerihelionTime());

        queryResults = searchTree.hyperRectangleSearch(queryPt, tolerances, myGeos);

        if (DEBUG) {
            std::cerr << "queryResults size: " << queryResults.size() << "  iteration: " << i << std::endl;
        }

        for(unsigned int j=0; j<queryResults.size(); j++){
            results.push_back(std::make_pair(queryResults.at(j).getValue(), i));
        }
    
        queryPt.clear();
    }

    return results;
}





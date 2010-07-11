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
 
/* File: orbitProximity.cc
 * Author: Matthew Cleveland
 * Purpose: 
 */

#include "lsst/mops/KDTree.h"
#include "lsst/mops/daymops/orbitProximity/orbitProximity.h"


#define DEBUG false


namespace lsst {
    namespace mops { 


// "internal" declarations

KDTree<unsigned int> buildKDTree(const std::vector<Orbit>);



/* use const by reference to avoid copying lots of memory */
std::vector<std::pair <unsigned int, unsigned int>  > 
getResults(const KDTree<unsigned int>&,
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
    KDTree<unsigned int> dataTree;
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
KDTree<unsigned int> buildKDTree(const std::vector<Orbit> orbits)
{
    KDTree<unsigned int> toReturn;
    
    std::vector<PointAndValue<unsigned int> > vecPV;
    
    const int orbitDimensions = 6;
    
    for(unsigned int i=0; i < orbits.size(); i++){
        
	PointAndValue<unsigned int> tempPV;
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
getResults(const KDTree<unsigned int> &searchTree,
           const std::vector<Orbit> &queryPoints,
           double maxPerihelion,
           double maxEccentricity,
           double maxInclination,
           double maxPerihelionArg,
           double maxLongitude,
           double maxPerihelionTime)
{
    std::vector<std::pair<unsigned int, unsigned int> > results;
    std::vector<PointAndValue<unsigned int> > queryResults;
    std::vector<double> queryPt;

    std::vector<GeometryType> myGeos;

    myGeos.push_back(CIRCULAR_DEGREES); //perihelion
    myGeos.push_back(EUCLIDEAN);        //eccentricity 
    myGeos.push_back(CIRCULAR_DEGREES); //inclination
    myGeos.push_back(CIRCULAR_DEGREES); //arg. of perihelion
    myGeos.push_back(CIRCULAR_DEGREES); //longitude
    myGeos.push_back(EUCLIDEAN);        //time of perihelion

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
            convertToStandardDegrees(queryPoints.at(i).getPerihelion()));
        //euclidean - eccentricity
        queryPt.push_back(queryPoints.at(i).getEccentricity());
        // degrees - inclination
        queryPt.push_back(
            convertToStandardDegrees(queryPoints.at(i).getInclination()));
        //degrees - arg. of perihelion (periapsis)
        queryPt.push_back(
            convertToStandardDegrees(queryPoints.at(i).getPerihelionArg()));
        //degrees - longitude
        queryPt.push_back(
            convertToStandardDegrees(queryPoints.at(i).getLongitude()));
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




    }} // close lsst::mops

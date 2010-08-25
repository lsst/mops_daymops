#define BOOST_TEST_MODULE orbitProximityTests

#include <boost/test/included/unit_test.hpp>
#include <boost/current_function.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>


#include "lsst/mops/Exceptions.h"
#include "lsst/mops/PointAndValue.h"
#include "lsst/mops/common.h"
#include "lsst/mops/KDTree.h"
#include "lsst/mops/Orbit.h"
#include "lsst/mops/daymops/orbitProximity/orbitProximity.h"


namespace lsst {
     namespace mops {

bool Eq(double a, double b) 
{
    double epsilon = 1e-10;
    return (fabs(a - b) < epsilon);
}



bool containsPair(unsigned int a, unsigned int b, std::vector<std::pair <unsigned int, unsigned int> > pairs) 
{     
  for (unsigned int i = 0; i < pairs.size(); i++) {
    
       std::pair <unsigned int, unsigned int> pair = pairs.at(i);
       
       if ((pair.first == a) && (pair.second == b)) {
	 return true;
       }
  }
  return false;
}






BOOST_AUTO_TEST_CASE( orbitProximity1 )
{
     
     // be a jerk, don't give any input
     std::vector<Orbit> dataOrbits;
     std::vector<Orbit> queryOrbits;
     std::vector<std::pair<unsigned int, unsigned int > > results = 
	  orbitProximity(dataOrbits, queryOrbits,
			 1, 1, 1, 1, 1, 1);
     BOOST_CHECK(results.size() == 0);
}


BOOST_AUTO_TEST_CASE( orbitProximity2 )
{
     
     // be a jerk, don't give any data
     std::vector<Orbit> dataOrbits;
     std::vector<Orbit> queryOrbits;
     Orbit tmpOrbit;
     tmpOrbit.setPerihelion(100);
     tmpOrbit.setEccentricity(100);
     tmpOrbit.setInclination(200);
     tmpOrbit.setPerihelionTime(53736);
     tmpOrbit.setLongitude(10);
     tmpOrbit.setEquinox(100);
     tmpOrbit.setOrbitID(1);
     queryOrbits.push_back(tmpOrbit);
     std::vector<std::pair<unsigned int, unsigned int > > results = 
	  orbitProximity(dataOrbits, queryOrbits,
			 1, 1, 1, 1, 1, 1);
     BOOST_CHECK(results.size() == 0);


}



BOOST_AUTO_TEST_CASE( orbitProximity3 )
{
     
     // be a jerk, don't give any query
     std::vector<Orbit> dataOrbits;
     std::vector<Orbit> queryOrbits;
     Orbit tmpOrbit;
     tmpOrbit.setPerihelion(100);
     tmpOrbit.setEccentricity(100);
     tmpOrbit.setInclination(200);
     tmpOrbit.setPerihelionTime(53736);
     tmpOrbit.setLongitude(10);
     tmpOrbit.setEquinox(100);
     tmpOrbit.setOrbitID(1);
     dataOrbits.push_back(tmpOrbit);
     std::vector<std::pair<unsigned int, unsigned int > > results = 
	  orbitProximity(dataOrbits, queryOrbits,
			 1, 1, 1, 1, 1, 1);
     BOOST_CHECK(results.size() == 0);
}




BOOST_AUTO_TEST_CASE( orbitProximity4 )
{
     
     // simple query - should match exactly one
     std::vector<Orbit> dataOrbits;
     std::vector<Orbit> queryOrbits;
     Orbit tmpOrbit;

     tmpOrbit.setPerihelion(100);
     tmpOrbit.setEccentricity(100);
     tmpOrbit.setInclination(200);
     tmpOrbit.setPerihelionTime(53736);
     tmpOrbit.setLongitude(10);
     tmpOrbit.setEquinox(100);
     tmpOrbit.setOrbitID(1);
     dataOrbits.push_back(tmpOrbit);

     tmpOrbit.setPerihelion(101);
     tmpOrbit.setEccentricity(101);
     tmpOrbit.setInclination(201);
     tmpOrbit.setPerihelionTime(53737);
     tmpOrbit.setLongitude(11);
     tmpOrbit.setEquinox(101);
     tmpOrbit.setOrbitID(2);
     queryOrbits.push_back(tmpOrbit);


     std::vector<std::pair<unsigned int, unsigned int > > results = 
	  orbitProximity(dataOrbits, queryOrbits,
			 1.1, 1.1, 1.1, 1.1, 1.1, 1.1);

     BOOST_CHECK(results.size() == 1);
     BOOST_CHECK(containsPair(0,0,results));
}



BOOST_AUTO_TEST_CASE( orbitProximity5 )
{
     
     // simple - should match exactly one
     std::vector<Orbit> dataOrbits;
     std::vector<Orbit> queryOrbits;
     Orbit tmpOrbit;

     tmpOrbit.setPerihelion(100);
     tmpOrbit.setEccentricity(100);
     tmpOrbit.setInclination(200);
     tmpOrbit.setPerihelionTime(53736);
     tmpOrbit.setLongitude(10);
     tmpOrbit.setEquinox(100);
     tmpOrbit.setOrbitID(1);
     dataOrbits.push_back(tmpOrbit);

     tmpOrbit.setPerihelion(101);
     tmpOrbit.setEccentricity(101);
     tmpOrbit.setInclination(201);
     tmpOrbit.setPerihelionTime(53737);
     tmpOrbit.setLongitude(11);
     tmpOrbit.setEquinox(101);
     tmpOrbit.setOrbitID(2);
     queryOrbits.push_back(tmpOrbit);

     tmpOrbit.setPerihelion(102); //perihelion too far!
     tmpOrbit.setEccentricity(101);
     tmpOrbit.setInclination(201);
     tmpOrbit.setPerihelionTime(53737);
     tmpOrbit.setLongitude(11);
     tmpOrbit.setEquinox(101);
     tmpOrbit.setOrbitID(2);
     queryOrbits.push_back(tmpOrbit);


     std::vector<std::pair<unsigned int, unsigned int > > results = 
	  orbitProximity(dataOrbits, queryOrbits,
			 1.1, 1.1, 1.1, 1.1, 1.1, 1.1);

     BOOST_CHECK(results.size() == 1);
     BOOST_CHECK(containsPair(0,0,results));


}







BOOST_AUTO_TEST_CASE( orbitProximity6 )
{
     
     // make sure perihelion, inclination, longitude, equinox are
     // all handled as degrees but perihelion time, eccentricity are not
     // - should match exactly one
     std::vector<Orbit> dataOrbits;
     std::vector<Orbit> queryOrbits;
     Orbit tmpOrbit;

     tmpOrbit.setPerihelion(0);
     tmpOrbit.setEccentricity(0);
     tmpOrbit.setInclination(0);
     tmpOrbit.setPerihelionTime(0);
     tmpOrbit.setLongitude(0);
     tmpOrbit.setEquinox(0);
     tmpOrbit.setOrbitID(1);
     dataOrbits.push_back(tmpOrbit);


     // this one should NOT match
     tmpOrbit.setPerihelion(359.5);
     tmpOrbit.setEccentricity(359.5); // eccentricity is NOT an angle
     tmpOrbit.setInclination(359.5);
     tmpOrbit.setPerihelionTime(0);
     tmpOrbit.setLongitude(359.5);
     tmpOrbit.setEquinox(359.5);
     tmpOrbit.setOrbitID(2);
     queryOrbits.push_back(tmpOrbit);


     // this one should NOT match
     tmpOrbit.setPerihelion(359.5);
     tmpOrbit.setEccentricity(0); 
     tmpOrbit.setInclination(359.5);
     tmpOrbit.setPerihelionTime(359.5); // perihelion time is NOT an angle
     tmpOrbit.setLongitude(359.5);
     tmpOrbit.setEquinox(359.5);
     tmpOrbit.setOrbitID(3);
     queryOrbits.push_back(tmpOrbit);

     // this one should match
     tmpOrbit.setPerihelion(359.5);
     tmpOrbit.setEccentricity(0);
     tmpOrbit.setInclination(359.5);
     tmpOrbit.setPerihelionTime(0);
     tmpOrbit.setLongitude(359.5);
     tmpOrbit.setEquinox(359.5);
     tmpOrbit.setOrbitID(4);
     queryOrbits.push_back(tmpOrbit);


     std::vector<std::pair<unsigned int, unsigned int > > results = 
	  orbitProximity(dataOrbits, queryOrbits,
				 1.1, 1.1, 1.1, 1.1, 1.1, 1.1);     

     BOOST_CHECK(results.size() == 1);
     BOOST_CHECK(containsPair(0,2,results));
     
}



     }} // close lsst::mops

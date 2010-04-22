// -*- LSST-C++ -*-
#define BOOST_TEST_MODULE detectionProximityUnitTests

#include <boost/test/included/unit_test.hpp>
#include <boost/current_function.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>


#include "lsst/mops/Exceptions.h"
#include "lsst/mops/daymops/detectionProximity/detectionProximity.h"


using namespace lsst::mops;


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






BOOST_AUTO_TEST_CASE( detectionProximity1 ) 
{
  // really simple example...
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  //           ID    MJD      RA      DEC     obs  NAME      MAG 
  Detection qd1(0, 53736,    100.00,   50.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd1(1, 53736,    100.01,   50.01, 566, "dummy1", 20.0, 0, 0);
  queryDets.push_back(qd1);
  dataDets.push_back(dd1);

  // these two are .0118 apart in angular distance
  queryResult = detectionProximity(queryDets, dataDets, .1, 1, 0.1);
  BOOST_CHECK(queryResult.size() == 1);
  BOOST_CHECK(containsPair(0,0,queryResult));

}




BOOST_AUTO_TEST_CASE( detectionProximity2 ) 
{
  // mean example - no data points
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  //           ID    MJD      RA      DEC     obs  NAME      MAG 
  Detection qd1(0, 53736,    100.00,   50.00, 566, "dummy1", 20.0, 0, 0);
  queryDets.push_back(qd1);

  queryResult = detectionProximity(queryDets, dataDets, .1, 1, 0.1);
  BOOST_CHECK(queryResult.size() == 0);

}




BOOST_AUTO_TEST_CASE( detectionProximity3 ) 
{
  // mean example - no query points
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  //           ID    MJD      RA      DEC     obs  NAME      MAG 
  Detection dd1(0, 53736,    100.00,   50.00, 566, "dummy1", 20.0, 0, 0);
  dataDets.push_back(dd1);

  queryResult = detectionProximity(queryDets, dataDets, .1, 1, 0.1);
  BOOST_CHECK(queryResult.size() == 0);

}





BOOST_AUTO_TEST_CASE( detectionProximity4 ) 
{
  // mean example - no points at all!
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  queryResult = detectionProximity(queryDets, dataDets, .1, 1, 0.1);
  BOOST_CHECK(queryResult.size() == 0);

}






BOOST_AUTO_TEST_CASE( detectionProximity5 ) 
{
  // pole crosser
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  //           ID    MJD      RA      DEC     obs  NAME      MAG 
  Detection qd1(0, 53736,    100.00,   89.50, 566, "dummy1", 20.0, 0, 0);
  Detection dd1(1, 53736,    280.00,   89.50, 566, "dummy1", 20.0, 0, 0);
  queryDets.push_back(qd1);
  dataDets.push_back(dd1);

  // these two are exactly 1.0 apart in angular distance
  queryResult = detectionProximity(queryDets, dataDets, 1.5, 1, 0.1);
  BOOST_CHECK(queryResult.size() == 1);
  BOOST_CHECK(containsPair(0,0,queryResult));

}




BOOST_AUTO_TEST_CASE( detectionProximity6 ) 
{
  // pole crosser, other pole
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  //           ID    MJD      RA      DEC     obs  NAME      MAG 
  Detection qd1(0, 53736,    100.00,   -89.50, 566, "dummy1", 20.0, 0, 0);
  Detection dd1(1, 53736,    280.00,   -89.50, 566, "dummy1", 20.0, 0, 0);
  queryDets.push_back(qd1);
  dataDets.push_back(dd1);

  // these two are exactly 1.0 apart in angular distance
  queryResult = detectionProximity(queryDets, dataDets, 1.5, 1, 0.1);
  BOOST_CHECK(queryResult.size() == 1);
  BOOST_CHECK(containsPair(0,0,queryResult));

}





BOOST_AUTO_TEST_CASE( detectionProximity7 ) 
{
  // simple, but with RA specified awkwardly
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  //           ID    MJD      RA      DEC     obs  NAME      MAG 
  Detection qd1(0, 53736,   -100.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd1(1, 53736,    260.01,   -10.00, 566, "dummy1", 20.0, 0, 0);
  queryDets.push_back(qd1);
  dataDets.push_back(dd1);

  // these two are ~.0098 apart in angular distance
  queryResult = detectionProximity(queryDets, dataDets, .1, 1, 0.1);
  BOOST_CHECK(queryResult.size() == 1);
  BOOST_CHECK(containsPair(0,0,queryResult));

}






BOOST_AUTO_TEST_CASE( detectionProximity8 ) 
{
  // lots of data points
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  //           ID    MJD      RA      DEC     obs  NAME      MAG 
  Detection qd1(0, 53736,    010.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd1(1, 53736,    007.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd2(2, 53736,    008.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd3(3, 53736,    009.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd4(4, 53736,    011.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd5(5, 53736,    012.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd6(6, 53736,    013.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  queryDets.push_back(qd1);
  dataDets.push_back(dd1);
  dataDets.push_back(dd2);
  dataDets.push_back(dd3);
  dataDets.push_back(dd4);
  dataDets.push_back(dd5);
  dataDets.push_back(dd6);

  // each is ~.98 apart in angular distance
  queryResult = detectionProximity(queryDets, dataDets, 1.01, 1, 0.1);
  BOOST_CHECK(queryResult.size() == 2);
  BOOST_CHECK(containsPair(0,2,queryResult));
  BOOST_CHECK(containsPair(0,3,queryResult));

}




BOOST_AUTO_TEST_CASE( detectionProximity9 ) 
{
  // lots of data points - and this time we test whether proper RA,Dec 
  // searching is being used or just euclidean approximation.
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  //           ID    MJD      RA      DEC     obs  NAME      MAG 
  Detection qd1(0, 53736,    010.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd1(1, 53736,    007.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd2(2, 53736,    008.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd3(3, 53736,    009.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd4(4, 53736,    011.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd5(5, 53736,    012.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd6(6, 53736,    013.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  queryDets.push_back(qd1);
  dataDets.push_back(dd1);
  dataDets.push_back(dd2);
  dataDets.push_back(dd3);
  dataDets.push_back(dd4);
  dataDets.push_back(dd5);
  dataDets.push_back(dd6);

  // each is ~.98 apart in angular distance
  queryResult = detectionProximity(queryDets, dataDets, .999, 1, 0.1);
  // if they used euclidean approximation, FAIL.
  BOOST_CHECK(queryResult.size() == 2);
  BOOST_CHECK(containsPair(0,2,queryResult));
  BOOST_CHECK(containsPair(0,3,queryResult));

}

BOOST_AUTO_TEST_CASE( detectionProximity10 ) 
{
  // lots of *query* points but just one data point
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  //           ID    MJD      RA      DEC     obs  NAME      MAG 
  Detection dd1(0, 53736,    010.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd1(1, 53736,    007.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd2(2, 53736,    008.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd3(3, 53736,    009.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd4(4, 53736,    011.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd5(5, 53736,    012.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd6(6, 53736,    013.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  dataDets.push_back(dd1);
  queryDets.push_back(qd1);
  queryDets.push_back(qd2);
  queryDets.push_back(qd3);
  queryDets.push_back(qd4);
  queryDets.push_back(qd5);
  queryDets.push_back(qd6);

  // each is ~.98 apart in angular distance
  queryResult = detectionProximity(queryDets, dataDets, 1.1, 1, 0.1);
  BOOST_CHECK(queryResult.size() == 2);
  BOOST_CHECK(containsPair(2,0,queryResult));
  BOOST_CHECK(containsPair(3,0,queryResult));

}


BOOST_AUTO_TEST_CASE( detectionProximity11 ) 
{
  // lots of *query* points but just one data point - again, test for
  // RA,Dec rather than euclidean approximation
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  //           ID    MJD      RA      DEC     obs  NAME      MAG 
  Detection dd1(0, 53736,    010.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd1(1, 53736,    007.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd2(2, 53736,    008.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd3(3, 53736,    009.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd4(4, 53736,    011.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd5(5, 53736,    012.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd6(6, 53736,    013.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  dataDets.push_back(dd1);
  queryDets.push_back(qd1);
  queryDets.push_back(qd2);
  queryDets.push_back(qd3);
  queryDets.push_back(qd4);
  queryDets.push_back(qd5);
  queryDets.push_back(qd6);

  // each is ~.98 apart in angular distance
  queryResult = detectionProximity(queryDets, dataDets, .999, 1, 0.1);
  // if they used euclidean approximation, FAIL.
  BOOST_CHECK(queryResult.size() == 2);
  BOOST_CHECK(containsPair(2,0,queryResult));
  BOOST_CHECK(containsPair(3,0,queryResult));

}





BOOST_AUTO_TEST_CASE( detectionProximity12 ) 
{
  // check that mag separation is honored
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  //           ID    MJD      RA      DEC     obs  NAME      MAG 
  Detection dd1(0, 53736,    010.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd1(1, 53736,    010.00,   -10.00, 566, "dummy1", 20.5, 0, 0);
  Detection qd2(1, 53736,    010.00,   -10.00, 566, "dummy1", 21.1, 0, 0);
  dataDets.push_back(dd1);
  queryDets.push_back(qd1);
  queryDets.push_back(qd2);

  queryResult = detectionProximity(queryDets, dataDets, .999, 1, 0.1);
  BOOST_CHECK(queryResult.size() == 1);
  BOOST_CHECK(containsPair(0,0,queryResult));
  
}





BOOST_AUTO_TEST_CASE( detectionProximity13 ) 
{
  // check that time separation is honored
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  //           ID    MJD      RA      DEC     obs  NAME      MAG 
  Detection dd1(0, 53736,   010.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd1(1, 53736.05,010.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection qd2(1, 53736.11,010.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  dataDets.push_back(dd1);
  queryDets.push_back(qd1);
  queryDets.push_back(qd2);

  queryResult = detectionProximity(queryDets, dataDets, .999, 1, 0.1);
  BOOST_CHECK(queryResult.size() == 1);
  BOOST_CHECK(containsPair(0,0,queryResult));
  
}





BOOST_AUTO_TEST_CASE( detectionProximity14 ) 
{
  // check that we catch RA 0/360 crossers
  std::vector<std::pair <unsigned int, unsigned int> > queryResult;
  std::vector<Detection> dataDets;
  std::vector<Detection> queryDets;
  
  /*Detection(long int ID, double epochMJD, double RA, double Dec, 
	    int obsCode, std::string objName, double mag,
	    double elongationLength, double elongationAngle); */

  //           ID    MJD      RA      DEC     obs  NAME      MAG 
  Detection qd1(0, 53736,   000.00,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd1(1, 53736,   359.99,   -10.00, 566, "dummy1", 20.0, 0, 0);
  Detection dd2(2, 53736,   360.01,   -10.00, 566, "dummy1", 20.0, 0, 0);
  queryDets.push_back(qd1);
  dataDets.push_back(dd1);
  dataDets.push_back(dd2);

  queryResult = detectionProximity(queryDets, dataDets, .999, 1, 0.1);
  BOOST_CHECK(queryResult.size() == 2);
  BOOST_CHECK(containsPair(0,0,queryResult));
  BOOST_CHECK(containsPair(0,1,queryResult));
  
}

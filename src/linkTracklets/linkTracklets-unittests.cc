// -*- LSST-C++ -*-
/* jonathan myers */
#define BOOST_TEST_MODULE linkTracklets

#include <boost/test/included/unit_test.hpp>
#include <boost/current_function.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>



#include "../Detection.h"
#include "../Tracklet.h"
#include "linkTracklets.h"

bool Eq(double a, double b) 
{
    double epsilon = 1e-10;
    return (fabs(a - b) < epsilon);
}






BOOST_AUTO_TEST_CASE( linkTracklets_whitebox_getBestFitVelocityAndAcceleration_test0) 
{
    std::vector<double> positions;
    std::vector<double> times;
    positions.push_back(0);
    positions.push_back(2);
    positions.push_back(6);
    times.push_back(0);
    times.push_back(1);
    times.push_back(2);
    double velocity, acceleration, position0;
    getBestFitVelocityAndAcceleration(positions, times, velocity, acceleration, position0);
    //std::cout << "position = " << position0 << " + " << velocity << "*t + " << acceleration << "*t^2" << std::endl;
    BOOST_CHECK(Eq(position0,    0));
    BOOST_CHECK(Eq(velocity,     1));
    BOOST_CHECK(Eq(acceleration, 1));    
}





BOOST_AUTO_TEST_CASE( linkTracklets_whitebox_getBestFitVelocityAndAcceleration_test1) 
{
    std::vector<double> positions;
    std::vector<double> times;
    positions.push_back(2);
    positions.push_back(6);
    positions.push_back(12);
    times.push_back(1);
    times.push_back(2);
    times.push_back(3);
    double velocity, acceleration, position0;
    getBestFitVelocityAndAcceleration(positions, times, velocity, acceleration, position0);
    //std::cout << "position = " << position0 << " + " << velocity << "*t + " << acceleration << "*t^2" << std::endl;
    BOOST_CHECK(Eq(position0,    0));
    BOOST_CHECK(Eq(velocity,     1));
    BOOST_CHECK(Eq(acceleration, 1));    
}






BOOST_AUTO_TEST_CASE( linkTracklets_blackbox_1 )
{
  // call with empty dets
  std::vector<Detection> myDets;
  std::vector<Tracklet> pairs;
  linkTrackletsConfig myConfig;
  std::vector<Track> results = linkTracklets(myDets, pairs, myConfig);
  BOOST_CHECK(pairs.size() == 0);
}




// helper function for creating sets of detections
void addDetectionAt(double MJD, double RA, double dec,  std::vector<Detection> &detVec)
{
    Detection tmpDet(detVec.size(), MJD, RA, dec, "566", "dummy",
                     24.0, 0., 0.);
    detVec.push_back(tmpDet);
}


void addPair(unsigned int id1, unsigned int id2, std::vector<Tracklet> &trackletVec) 
{
    Tracklet tmpTracklet;
    tmpTracklet.indices.insert(id1);
    tmpTracklet.indices.insert(id2);
    trackletVec.push_back(tmpTracklet);
}



BOOST_AUTO_TEST_CASE( linkTracklets_easy_1 )
{
  std::vector<Detection> myDets;
  addDetectionAt(5300.0,  50,     50, myDets);
  addDetectionAt(5300.01, 50.001, 50.001, myDets);
  addDetectionAt(5301.0,  50.1,   50.1, myDets);
  addDetectionAt(5301.01, 50.101, 50.101, myDets);
  addDetectionAt(5302.0,  50.2,   50.2, myDets);
  addDetectionAt(5302.01, 50.201, 50.201, myDets);


  std::vector<Tracklet> pairs;
  addPair(0,1, pairs);
  addPair(2,3, pairs);
  addPair(4,5, pairs);
  
  linkTrackletsConfig myConfig;

  std::vector<Track> results = linkTracklets(myDets, pairs, myConfig);

  BOOST_CHECK(results.size() == 1);
}





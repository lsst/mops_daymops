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






// helper function for creating sets of detections
void addDetectionAt(double MJD, double RA, double dec,  std::vector<Detection> &detVec)
{
    Detection tmpDet(detVec.size(), MJD, RA, dec, 566, "dummy",
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


BOOST_AUTO_TEST_CASE( linkTracklets_blackbox_1 )
{
  // call with empty dets
  std::vector<Detection> myDets;
  std::vector<Tracklet> pairs;
  linkTrackletsConfig myConfig;
  std::vector<Track> results = linkTracklets(myDets, pairs, myConfig);
  BOOST_CHECK(pairs.size() == 0);
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







BOOST_AUTO_TEST_CASE( linkTracklets_easy_2 )
{
    // same as 1, but with more tracks (all clearly separated)

  std::vector<Detection> myDets;
  std::vector<Tracklet> pairs;
  for (unsigned int i = 0; i < 10; i++) {

      addDetectionAt(5300.0,  50 + i,     50, myDets);
      addDetectionAt(5300.01, 50.001 + i, 50.001, myDets);
      addDetectionAt(5301.0,  50.1 + i,   50.1, myDets);
      addDetectionAt(5301.01, 50.101 + i, 50.101, myDets);
      addDetectionAt(5302.0,  50.2 + i,   50.2, myDets);
      addDetectionAt(5302.01, 50.201 + i, 50.201, myDets);

      addPair(0 + 6*i,1 + 6*i, pairs);
      addPair(2 + 6*i,3 + 6*i, pairs);
      addPair(4 + 6*i,5 + 6*i, pairs);
      
  }

  
  linkTrackletsConfig myConfig;

  std::vector<Track> results = linkTracklets(myDets, pairs, myConfig);

  BOOST_CHECK(results.size() == 10);
}






BOOST_AUTO_TEST_CASE( linkTracklets_easy_3 )
{
    // same as 2, but with tracks crossing RA 0 line

  std::vector<Detection> myDets;
  std::vector<Tracklet> pairs;
  for (unsigned int i = 0; i < 10; i++) {

      addDetectionAt(5300.0,  50 + i,     50, myDets);
      addDetectionAt(5300.01, 50.001 + i, 50.001, myDets);
      addDetectionAt(5301.0,  50.1 + i,   50.1, myDets);
      addDetectionAt(5301.01, 50.101 + i, 50.101, myDets);
      addDetectionAt(5302.0,  50.2 + i,   50.2, myDets);
      addDetectionAt(5302.01, 50.201 + i, 50.201, myDets);

      addPair(0 + 6*i,1 + 6*i, pairs);
      addPair(2 + 6*i,3 + 6*i, pairs);
      addPair(4 + 6*i,5 + 6*i, pairs);
      
  }

  
  linkTrackletsConfig myConfig;

  std::vector<Track> results = linkTracklets(myDets, pairs, myConfig);

  BOOST_CHECK(results.size() == 10);
}





BOOST_AUTO_TEST_CASE( linkTracklets_easy_4_1 )
{
    // same as 1, but with track crossing RA 0 line

  std::vector<Detection> myDets;
  std::vector<Tracklet> pairs;

  addDetectionAt(5300.0,  359.9,       50, myDets);
  addDetectionAt(5300.01, 359.901, 50.001, myDets);
  addDetectionAt(5301.0,  0.,        50.1, myDets);
  addDetectionAt(5301.01, 0.001,    50.101, myDets);
  addDetectionAt(5302.0,   0.1,      50.2, myDets);
  addDetectionAt(5302.01,  0.101,  50.201, myDets);
  
  addPair(0,1, pairs);
  addPair(2,3, pairs);
  addPair(4,5, pairs);
  
  linkTrackletsConfig myConfig;

  std::vector<Track> results = linkTracklets(myDets, pairs, myConfig);

  //std::cout << "results size = " << results.size() << '\n';
  BOOST_CHECK(results.size() == 1);
}



BOOST_AUTO_TEST_CASE( linkTracklets_easy_4 )
{
    // same as 2, but with tracks crossing RA 0 line

  std::vector<Detection> myDets;
  std::vector<Tracklet> pairs;
  for (unsigned int i = 0; i < 10; i++) {

      addDetectionAt(5300.0,  359.9,       50 + i, myDets);
      addDetectionAt(5300.01, 359.901, 50.001 + i, myDets);
      addDetectionAt(5301.0,  0.,        50.1 + i, myDets);
      addDetectionAt(5301.01, 0.001,   50.101 + i, myDets);
      addDetectionAt(5302.0,   0.1,      50.2 + i, myDets);
      addDetectionAt(5302.01,  0.101,  50.201 + i, myDets);

      addPair(0 + 6*i,1 + 6*i, pairs);
      addPair(2 + 6*i,3 + 6*i, pairs);
      addPair(4 + 6*i,5 + 6*i, pairs);
      
  }

  
  linkTrackletsConfig myConfig;

  std::vector<Track> results = linkTracklets(myDets, pairs, myConfig);

  //std::cout << "results size = " << results.size() << '\n';
  BOOST_CHECK(results.size() == 10);
}



BOOST_AUTO_TEST_CASE( linkTracklets_easy_5 )
{
    // same as 4, but going the other way!

  std::vector<Detection> myDets;
  std::vector<Tracklet> pairs;
  for (unsigned int i = 0; i < 10; i++) {

      addDetectionAt(5302.0,  0.101,  50.201 + i, myDets);
      addDetectionAt(5302.01,   0.1,      50.2 + i, myDets);
      addDetectionAt(5301.0, 0.001,   50.101 + i, myDets);
      addDetectionAt(5301.01,  0.,        50.1 + i, myDets);
      addDetectionAt(5300.0, 359.901, 50.001 + i, myDets);
      addDetectionAt(5300.01,  359.9,       50 + i, myDets);

      addPair(0 + 6*i,1 + 6*i, pairs);
      addPair(2 + 6*i,3 + 6*i, pairs);
      addPair(4 + 6*i,5 + 6*i, pairs);
      
  }

  
  linkTrackletsConfig myConfig;

  std::vector<Track> results = linkTracklets(myDets, pairs, myConfig);

  //std::cout << "results size = " << results.size() << '\n';
  BOOST_CHECK(results.size() == 10);
}

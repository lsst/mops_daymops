// -*- LSST-C++ -*-
/* jonathan myers */
#define BOOST_TEST_MODULE findTracklets

#include <boost/test/included/unit_test.hpp>
#include <boost/current_function.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>



#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/Tracklet.h"
#include "lsst/mops/PointAndValue.h"
#include "lsst/mops/common.h"
#include "lsst/mops/KDTree.h"
#include "lsst/mops/daymops/findTracklets/findTracklets.h"


using namespace lsst::mops;

bool Eq(double a, double b) 
{
    double epsilon = 1e-10;
    return (fabs(a - b) < epsilon);
}

BOOST_AUTO_TEST_CASE( findTracklets_blackbox_1 )
{
  // call with empty dets
  std::vector<MopsDetection> myDets;
  TrackletVector *pairs;
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 1.5;
  pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 0);
  delete pairs;
}


// helper function for creating sets of detections
void addDetectionAt(double MJD, double RA, double dec,  std::vector<MopsDetection> &detVec)
{
    MopsDetection tmpDet(detVec.size(), MJD, RA, dec);
    detVec.push_back(tmpDet);
}


//helper function for grading correctness

bool containsPair(unsigned int ID1, unsigned int ID2, TrackletVector *pairs)
{
    for (unsigned int i = 0; i < pairs->size(); i++) {
        Tracklet curPair = pairs->at(i);
        if ((curPair.indices.find(ID1) != curPair.indices.end())
            && (curPair.indices.find(ID2) != curPair.indices.end()))
            return true;
    }
    return false;
}


// test that tracklet of 2 dets is created
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_2 )
{
  std::vector<MopsDetection> myDets;
  addDetectionAt(53736.0, 100.0, 10.0, myDets);
  addDetectionAt(53737.0, 100.1, 10.1, myDets);

  //std::cerr << "boost test 2" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 1.5;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  delete pairs;
}





// test 3 dets of same object gets 3 pairs
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_3 )
{
  std::vector<MopsDetection> myDets;
  addDetectionAt(53736.0, 100.0, 10.0, myDets); // id 0
  addDetectionAt(53737.0, 100.1, 10.1, myDets); // id 1
  addDetectionAt(53738.0, 100.2, 10.2, myDets); // id 2
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;
  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 3);
  BOOST_CHECK(containsPair(0, 1, pairs));
  BOOST_CHECK(containsPair(0, 2, pairs));
  BOOST_CHECK(containsPair(1, 2, pairs));
  delete pairs;
}



// 2 image times, 2 objects = 2 tracklets
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_4 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 100.0, -10.0, myDets); // id 0
  addDetectionAt(53737.0, 100.1, -10.1, myDets); // id 1
  // second "object"
  addDetectionAt(53736.0, 150.0, -10.0, myDets); // id 2
  addDetectionAt(53737.0, 150.1, -10.1, myDets); // id 3
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;
  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 2);
  BOOST_CHECK(containsPair(0, 1, pairs));
  BOOST_CHECK(containsPair(2, 3, pairs));
  delete pairs;
}





// check that velocity filter is working
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_5 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 100.0, 80.0, myDets); // id 0
  addDetectionAt(53737.0, 100.1, 80.1, myDets); // id 1
  // second "object" - too fast! (1.104 deg/day)
  addDetectionAt(53736.0, 150.0, 0.0, myDets); // id 2
  addDetectionAt(53737.0, 151.1, 0.1, myDets); // id 3
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;
  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}


// check that velocity filter is working (for time separation != 1 day)
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_5_1 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 100.0, -80.0, myDets); // id 0
  addDetectionAt(53736.5, 100.02, -80.1, myDets); // id 1
  // second "object" - too fast! (1.104 deg/day)
  addDetectionAt(53736.0, 150.0, -0.0, myDets); // id 3
  addDetectionAt(53736.5, 150.55, -0.1, myDets); // id 4
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;
  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}

// check that velocity filter is working, with things moving too fast in dec
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_5_2 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 100.0, -80.0, myDets); // id 0
  addDetectionAt(53736.5, 100.02, -80.1, myDets); // id 1
  // second "object" - too fast! (1.1 deg/day)
  addDetectionAt(53736.0, 150.0, -0.0, myDets); // id 3
  addDetectionAt(53736.5, 150.0, -0.55, myDets); // id 4
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;
  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}

// check that velocity filter is working, with things moving too fast in RA near a pole
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_5_3 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 100.0, -80.0, myDets); // id 0
  addDetectionAt(53736.5, 100.02, -80.1, myDets); // id 1
  // second "object" - too fast! (~1.0418 deg/day)
  addDetectionAt(53736.0, 150.0, -80.0, myDets); // id 3
  addDetectionAt(53736.5, 153.0, -80.0, myDets); // id 4
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;
  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}




// check that velocity filter is working *correctly* - that we search in a circle, not a square.
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_6 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 100.0, 20.0, myDets); // id 0
  addDetectionAt(53737.0, 100.1, 20.1, myDets); // id 1
  // second "object" - too fast! ~1.23 deg/day
  addDetectionAt(53736.0, 150.0, 20.0, myDets); // id 3
  addDetectionAt(53737.0, 150.9, 20.9, myDets); // id 4

  //  std::cerr << "boost test 6" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;
  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}


// check that velocity filter is working *correctly* - that we search in a circle, not a square.
// and also that this works for time separation != 1 day
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_6_1 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 100.0, 60.0, myDets); // id 0
  addDetectionAt(53736.5, 100.05, 60.05, myDets); // id 1
  // second "object" - too fast! ~1.004 deg/day
  addDetectionAt(53736.0, 150.0, 60.0, myDets); // id 3
  addDetectionAt(53736.5, 150.45, 60.45, myDets); // id 4

  //  std::cerr << "boost test 6" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;
  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}




// see if we catch an object crossing the RA 0/360 line
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_7 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 359.9, 60.0, myDets); // id 0
  addDetectionAt(53737.0, 000.1, 60.1, myDets); // id 1
  //  std::cerr << "boost test 7" << std::endl;
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;
  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}


// see if we catch an object crossing the north pole
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_8 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 0.0, 89.9, myDets); // id 0
  addDetectionAt(53737.0, 180.0, 89.9, myDets); // id 1

  //  std::cerr << "boost test 8" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;
  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}



// see if we catch an object crossing the RA 0/360 line (backwards)
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_9 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 000.1, 50.0, myDets); // id 0
  addDetectionAt(53737.0, 359.9, 50.1, myDets); // id 1

  //  std::cerr << "boost test 9" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;
  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}


// see if we catch an object crossing the north pole
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_10 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 100.0, 89.9, myDets); // id 0
  addDetectionAt(53737.0, 280.0, 89.9, myDets); // id 1

  //std:cerr << "boost test 10" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;
  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}

// see if we catch an object crossing the south pole (backwards)
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_11 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 100.0, -89.9, myDets); // id 0
  addDetectionAt(53737.0, 280.0, -89.9, myDets); // id 1

  //std:cerr << "boost test 10" << std::endl;
 
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}


// non-standard dec values, crossing south pole!
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_112)
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 100.0, -89.9, myDets); // id 0
  addDetectionAt(53737.0, 280.0, -90.1, myDets); // id 1

  //std:cerr << "boost test 10" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}






// check that velocity is still handled correctly on pole crossing
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_13 )
{
  std::vector<MopsDetection> myDets;
  // first "object" - too fast for 1 deg/day velocity filter (1.1 deg/day)
  addDetectionAt(53736.0, 180.45, 89.45, myDets); // id 0
  addDetectionAt(53737.0, 000.45, 89.45, myDets); // id 1

  //std:cerr << "boost test 13" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 0);
  delete pairs;
}



// check that negative RA, Dec are interpreted correctly
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_14 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, -180.1, -0.1, myDets); // id 0
  addDetectionAt(53737.0, 180.0, 359.9, myDets); // id 1

  // second "object" -- too fast in RA (~1.05 deg/day)
  addDetectionAt(53736.0, -283.0, -280.9, myDets); // id 3
  addDetectionAt(53737.0, 80.0, 80.0, myDets); // id 4

  //std:cerr << "boost test 14" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}






// one really brain-dead test: just do lots and lots of objects
// lots and lots of noise
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_15 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 100.0, 100.0, myDets); // id 0
  addDetectionAt(53737.0, 100.5, 100.0, myDets); // id 1

  // three noise points in first 'image'
  addDetectionAt(53736.0, 10.0, 20.0, myDets); // id 2
  addDetectionAt(53736.0, 180.0, 70.0, myDets); // id 3
  addDetectionAt(53736.0, 200.0, 20.0, myDets); // id 4

  // second "object"
  addDetectionAt(53736.0, 5.0,  100.0, myDets); // id 5
  addDetectionAt(53737.0, 5.01, 99.99, myDets); // id 6

  // three noise points in second image
  addDetectionAt(53737.0, 12.0, 20.0, myDets); // id 7
  addDetectionAt(53737.0, 108.0, -20.0, myDets); // id 8
  addDetectionAt(53737.0, 50.0, -80.0, myDets); // id 9

  // third "object"
  addDetectionAt(53736.0, 355.0,  5.0, myDets); // id 10
  addDetectionAt(53737.0, 354.9, 5.1, myDets); // id 11

  // three more noise points in first image
  addDetectionAt(53736.0, 250.0, 30.0, myDets); // id 12
  addDetectionAt(53736.0, 280.0, -30.0, myDets); // id 13
  addDetectionAt(53736.0, 001.0, -40.0, myDets); // id 14

  // fourth "object"
  addDetectionAt(53736.0, 210.0,  -15.0, myDets); // id 15
  addDetectionAt(53737.0, 210.0, -15.5, myDets); // id 16

  // three more noise points in first image
  addDetectionAt(53737.0, 005.0, 005.0, myDets); // id 17
  addDetectionAt(53737.0, 010.0, 005.0, myDets); // id 18
  addDetectionAt(53737.0, 290.0, 20.0, myDets); // id 19


  //std:cerr << "boost test 15" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 4);
  BOOST_CHECK(containsPair(0,  1, pairs));
  BOOST_CHECK(containsPair(5,  6, pairs));
  BOOST_CHECK(containsPair(10, 11, pairs));
  BOOST_CHECK(containsPair(15, 16, pairs));
  delete pairs;

}


//make sure we don't link detections from the same image.
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_16 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 100.0, 100.0, myDets); // id 0
  addDetectionAt(53736.0, 100.1, 100.0, myDets); // id 1

  //std:cerr << "boost test 16" << std::endl;
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 0);
  delete pairs;
}


// test that negative RA is handled sanely
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_17 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 000.1, 0.0, myDets); // id 0
  addDetectionAt(53737.0, -00.1, 0.0, myDets); // id 1
  // noise
  addDetectionAt(53736.0, -01.5, 0.0, myDets); //id 2

  //std:cerr << "boost test 17" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}


// test that negative Dec is handled sanely
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_18 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 000.1, -00.1, myDets); // id 0
  addDetectionAt(53737.0, 000.1, 000.1, myDets); // id 1
  // noise
  addDetectionAt(53736.0, 000.1, -02.0, myDets); //id 2

  //std:cerr << "boost test 18" << std::endl;
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(0, 1, pairs));
  delete pairs;
}



// check that behavior is sane when no results are returned
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_19 )
{
  std::vector<MopsDetection> myDets;
  // only "object" - too fast!
  addDetectionAt(53736.0, 150.0, 0.0, myDets); // id 0
  addDetectionAt(53737.0, 151.1, 0.1, myDets); // id 1

  //std:cerr << "boost test 19" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 0);
  delete pairs;
}




// /*********************************************************************/
// /*********************************************************************/
// /********** tests with minv != 0 *************************************/
// /*********************************************************************/
// /*********************************************************************/


// test that minv does something in RA
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_using_minv_1 )
{
  std::vector<MopsDetection> myDets;
  // first "object" -- too slow
  addDetectionAt(53736.0, 100.0, 100.0, myDets); // id 0
  addDetectionAt(53737.0, 100.1, 100.0, myDets); // id 1
  
  // second "object", .689 deg/day
  addDetectionAt(53736.0, 300.0, 10.0, myDets); // id 2
  addDetectionAt(53737.0, 300.7, 10.0, myDets); // id 3

  //std:cerr << "boost min test 1" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.minV = .5;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(2, 3, pairs));
  delete pairs;
  
}


// test that minv does something in RA
// (again, with time separation != exactly 1 day)
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_using_minv_1_1 )
{
  std::vector<MopsDetection> myDets;
  // first "object" -- too slow (.017 deg/day)
  addDetectionAt(53736.0, 100.0, 80.0, myDets); // id 0
  addDetectionAt(53736.5, 100.05, 80.0, myDets); // id 1
  
  // second "object" (.59 deg/day)
  addDetectionAt(53736.0, 300.0, -10.0, myDets); // id 2
  addDetectionAt(53736.5, 300.3, -10.0, myDets); // id 3

  //std:cerr << "boost min test 1" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.minV = .5;
  config.maxDt = 3.0;
  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(2, 3, pairs));
  delete pairs;
  
}



// test that minv does something in dec
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_using_minv_2 )
{
  std::vector<MopsDetection> myDets;
  // first "object" -- too slow
  addDetectionAt(53736.0, 100.0, -10.0, myDets); // id 0
  addDetectionAt(53737.0, 100.0, -10.1, myDets); // id 1
  
  // second "object" 
  addDetectionAt(53736.0, 300.0, 80.0, myDets); // id 2
  addDetectionAt(53737.0, 300.0, 80.6, myDets); // id 3

  //std:cerr << "boost min test 2" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.minV = .5;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(2, 3, pairs));
  delete pairs;
  
}


// test that minv does something in dec
// again, with time separation != exactly one day
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_using_minv_2_1 )
{
  std::vector<MopsDetection> myDets;
  // first "object" -- too slow
  addDetectionAt(53736.0, 100.0, 100.0, myDets); // id 0
  addDetectionAt(53736.5, 100.0, 100.05, myDets); // id 1
  
  // second "object" 
  addDetectionAt(53736.0, 300.0, 300.0, myDets); // id 2
  addDetectionAt(53736.5, 300.0, 300.3, myDets); // id 3

  //std:cerr << "boost min test 2" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.minV = .5;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(2, 3, pairs));
  delete pairs;
  
}


// test that minv does something in RA.dec
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_using_minv_3 )
{
  std::vector<MopsDetection> myDets;
  // first "object" -- too slow (distance of ~.14)
  addDetectionAt(53736.0, 100.0, 10.0, myDets); // id 0
  addDetectionAt(53737.0, 100.1, 10.1, myDets); // id 1
  
  // second "object" -- just fast enough (distance of ~.66)
  addDetectionAt(53736.0, 300.0, 30.0, myDets); // id 2
  addDetectionAt(53737.0, 300.5, 30.5, myDets); // id 3

  //std:cerr << "boost min test 3" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.minV = .5;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(2, 3, pairs));
  delete pairs;
  
}




// test that minv works with crossing RA 0 line
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_using_minv_4 )
{
  std::vector<MopsDetection> myDets;
  // first "object" -- too slow
  addDetectionAt(53736.0, 000.05, 000.05, myDets); // id 0
  addDetectionAt(53737.0, 359.95, 000.05, myDets); // id 1
  
  // second "object" -- just fast enough
  addDetectionAt(53736.0, 000.1, 020.5, myDets); // id 2
  addDetectionAt(53737.0, 359.5, 020.5, myDets); // id 3

  //std:cerr << "boost min test 4" << std::endl;
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.minV = .5;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(2, 3, pairs));
  delete pairs;
  
}





// test that minv works with crossing north pole
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_using_minv_5 )
{
  std::vector<MopsDetection> myDets;
  // first "object" -- too slow (distance of .1)
  addDetectionAt(53736.0, 000.05, 89.95, myDets); // id 0
  addDetectionAt(53737.0, 180.05, 89.95, myDets); // id 1
  
  // second "object" -- just fast enough (distance of .6)
  addDetectionAt(53736.0, 010.5, 89.7, myDets); // id 2
  addDetectionAt(53737.0, 190.5, 89.7, myDets); // id 3

  //std:cerr << "boost min test 5" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.minV = .5;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(2, 3, pairs));
  delete pairs;
  
}





// test that minv works with crossing RA *and* dec 0 lines
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_using_minv_6 )
{
  std::vector<MopsDetection> myDets;
  // first "object" -- too slow (distance of ~.14)
  addDetectionAt(53736.0, 000.05, 000.05, myDets); // id 0
  addDetectionAt(53737.0, 359.95, -00.05, myDets); // id 1
  
  // second "object" -- just fast enough (distance of ~.707)
  addDetectionAt(53736.0, 000.25, 000.25, myDets); // id 2
  addDetectionAt(53737.0, 359.75, -00.25, myDets); // id 3

  //std:cerr << "boost min test 6" << std::endl;

  findTrackletsConfig config;
  config.maxV = 1.0;
  config.minV = .5;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 1);
  BOOST_CHECK(containsPair(2, 3, pairs));
  delete pairs;
  
}







// test that behavior is sane when no tracklets are returned
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_using_minv_7 )
{
  std::vector<MopsDetection> myDets;
  // first "object" -- too slow (distance of ~.14)
  addDetectionAt(53736.0, 000.05, 000.05, myDets); // id 0
  addDetectionAt(53737.0, 359.05, 359.05, myDets); // id 1
  
  // second "object" -- also too slow(distance of ~.707)
  addDetectionAt(53736.0, 000.25, 000.25, myDets); // id 2
  addDetectionAt(53737.0, 359.75, 359.75, myDets); // id 3

  //std:cerr << "boost min test 7" << std::endl;
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.minV = .8;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 0);
  delete pairs;

}







// one more really brain-dead test: just do lots and lots of objects
// lots and lots of noise
BOOST_AUTO_TEST_CASE( findTracklets_blackbox_using_minv_8 )
{
  std::vector<MopsDetection> myDets;
  // first "object"
  addDetectionAt(53736.0, 100.0, 10.0, myDets); // id 0
  addDetectionAt(53737.0, 100.5, 10.0, myDets); // id 1

  // three noise points in first image
  addDetectionAt(53736.0, 10.0, 20.0, myDets); // id 2
  addDetectionAt(53736.0, 180.0, 100.0, myDets); // id 3
  addDetectionAt(53736.0, 200.0, 20.0, myDets); // id 4

  // second "object" -- too slow!
  addDetectionAt(53736.0, 5.0,  10.0, myDets); // id 5
  addDetectionAt(53737.0, 5.01,  9.99, myDets); // id 6

  // three noise points in second image
  addDetectionAt(53737.0, 12.0, 20.0, myDets); // id 7
  addDetectionAt(53737.0, 108.0, -20.0, myDets); // id 8
  addDetectionAt(53737.0, 50.0, 350.0, myDets); // id 9

  // third "object" - too slow, only ~.14 apart
  addDetectionAt(53736.0, 355.0,  5.0, myDets); // id 10
  addDetectionAt(53737.0, 354.9, 5.1, myDets); // id 11

  // three more noise points in first image
  addDetectionAt(53736.0, 250.0, 80.0, myDets); // id 12
  addDetectionAt(53736.0, 280.0, -80.0, myDets); // id 13
  addDetectionAt(53736.0, 001.0, 20.0, myDets); // id 14

  // fourth "object"
  addDetectionAt(53736.0, 210.0,  15.0, myDets); // id 15
  addDetectionAt(53737.0, 210.0, 15.5, myDets); // id 16

  // three more noise points in first image
  addDetectionAt(53737.0, 005.0, 005.0, myDets); // id 17
  addDetectionAt(53737.0, 010.0, 005.0, myDets); // id 18
  addDetectionAt(53737.0, 290.0, 20.0, myDets); // id 19

  //std:cerr << "boost min test 8" << std::endl;
  findTrackletsConfig config;
  config.maxV = 1.0;
  config.minV = .2;
  config.maxDt = 3.0;

  TrackletVector *pairs = findTracklets(myDets, config);
  BOOST_CHECK(pairs->size() == 2);
  BOOST_CHECK(containsPair(0,  1, pairs));
  //BOOST_CHECK(containsPair(5,  6, pairs));
  //BOOST_CHECK(containsPair(10, 11, pairs));
  BOOST_CHECK(containsPair(15, 16, pairs));
  delete pairs;

}

// -*- LSST-C++ -*-

/* jmyers  6/1/2011
 *
 * the new fieldProximity demands two things: first, that tracks
 * contain 1 detection per day, and second, that the track contains at
 * least one detection within the square which enscribes the
 * (circular) image.  This will only break down if tracks move very
 * quickly - by my calculations, >1.44 deg/day given an image of
 * radius 1.75 degrees.
 *
 */ 

#define BOOST_TEST_MODULE fieldProximityTests

#include <boost/test/included/unit_test.hpp>
#include <boost/current_function.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>


#include "lsst/mops/Exceptions.h"
#include "lsst/mops/daymops/fieldProximity/Field.h"
#include "lsst/mops/daymops/fieldProximity/fieldProximity.h"
#include "lsst/mops/daymops/fieldProximity/TrackForFieldProximity.h"






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






BOOST_AUTO_TEST_CASE ( fieldProximity1 ) 
{
     //simple test - one field, one track.
     std::vector<Field> queryFields;
     std::vector<FieldProximityTrack> allTracks;

     Field tmpField;
     tmpField.setFieldID(1);
     tmpField.setEpochMJD(300);
     tmpField.setRA(50);
     tmpField.setDec(50);
     tmpField.setRadius(1.75);
     queryFields.push_back(tmpField);
     
     FieldProximityTrack tmpTrack;
     tmpTrack.setID(42);
     FieldProximityPoint tmpPoint;
     tmpPoint.setRA(50);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(299.5);
     tmpTrack.addPoint(tmpPoint);

     tmpPoint.setRA(51);
     tmpPoint.setDec(49);
     tmpPoint.setEpochMJD(300.5);
     tmpTrack.addPoint(tmpPoint);

     allTracks.push_back(tmpTrack);

     std::vector<std::pair <unsigned int, unsigned int> > pairs =
	  fieldProximity(allTracks, queryFields, 0);
     
     BOOST_CHECK(pairs.size() == 1);
     BOOST_CHECK(containsPair(1,42,pairs));
     
}




BOOST_AUTO_TEST_CASE ( fieldProximity2 ) 
{
     //simple test - one field, two tracks.
     std::vector<Field> queryFields;
     std::vector<FieldProximityTrack> allTracks;

     Field tmpField;
     tmpField.setFieldID(1);
     tmpField.setEpochMJD(300);
     tmpField.setRA(50);
     tmpField.setDec(50);
     tmpField.setRadius(1.75);
     queryFields.push_back(tmpField);
     
     FieldProximityTrack tmpTrack;
     tmpTrack.setID(42);
     FieldProximityPoint tmpPoint;
     tmpPoint.setRA(50);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(299.5);
     tmpTrack.addPoint(tmpPoint);

     tmpPoint.setRA(51);
     tmpPoint.setDec(49);
     tmpPoint.setEpochMJD(300.5);
     tmpTrack.addPoint(tmpPoint);
     allTracks.push_back(tmpTrack);
     
     //note that this one is too far west
     FieldProximityTrack tmpTrack2;
     tmpTrack2.setID(33);
     tmpPoint.setRA(51.76);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(299.5);
     tmpTrack2.addPoint(tmpPoint);

     tmpPoint.setRA(52.5);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(300.5);
     tmpTrack2.addPoint(tmpPoint);
     allTracks.push_back(tmpTrack2);     

     std::vector<std::pair <unsigned int, unsigned int> > pairs =
	  fieldProximity(allTracks, queryFields, 0);
     
     BOOST_CHECK(pairs.size() == 1);
     BOOST_CHECK(containsPair(1,42,pairs));
     
}






BOOST_AUTO_TEST_CASE ( fieldProximity4 ) 
{
     //try to cause trouble - give no tracks.
     std::vector<Field> queryFields;
     std::vector<FieldProximityTrack> allTracks;

     Field tmpField;
     tmpField.setFieldID(1);
     tmpField.setEpochMJD(300);
     tmpField.setRA(50);
     tmpField.setDec(50);
     tmpField.setRadius(10);
     queryFields.push_back(tmpField);
     
     std::vector<std::pair <unsigned int, unsigned int> > pairs =
	  fieldProximity(allTracks, queryFields, 0);
     
     BOOST_CHECK(pairs.size() == 0);
     //BOOST_CHECK(containsPair(0,0,pairs));
     
}



BOOST_AUTO_TEST_CASE ( fieldProximity5 ) 
{
     //try to cause trouble - give nothing at all.
     std::vector<Field> queryFields;
     std::vector<FieldProximityTrack> allTracks;     

     std::vector<std::pair <unsigned int, unsigned int> > pairs =
	  fieldProximity(allTracks, queryFields, 0);
     
     BOOST_CHECK(pairs.size() == 0);
     //BOOST_CHECK(containsPair(0,0,pairs));
     
}



BOOST_AUTO_TEST_CASE ( fieldProximity6 ) 
{
     //try to cause trouble - give no fields!
     std::vector<Field> queryFields;
     std::vector<FieldProximityTrack> allTracks;     

     FieldProximityTrack tmpTrack;
     FieldProximityPoint tmpPoint;
     tmpPoint.setRA(30);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(299.5);
     tmpTrack.addPoint(tmpPoint);

     tmpPoint.setRA(70);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(300.5);
     tmpTrack.addPoint(tmpPoint);
     allTracks.push_back(tmpTrack);

     std::vector<std::pair <unsigned int, unsigned int> > pairs =
	  fieldProximity(allTracks, queryFields, 0);
     
     BOOST_CHECK(pairs.size() == 0);
     //BOOST_CHECK(containsPair(0,0,pairs));
     
}





BOOST_AUTO_TEST_CASE ( fieldProximity7 ) 
{
     //tricky test - make the track cross the 0/360 line in RA..
     std::vector<Field> queryFields;
     std::vector<FieldProximityTrack> allTracks;

     Field tmpField;
     tmpField.setFieldID(1);
     tmpField.setEpochMJD(300);
     tmpField.setRA(1.5);
     tmpField.setDec(50);
     tmpField.setRadius(1.75);
     queryFields.push_back(tmpField);
     
     FieldProximityTrack tmpTrack;
     FieldProximityPoint tmpPoint;
     tmpTrack.setID(42);
     tmpPoint.setRA(359);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(299.5);
     tmpTrack.addPoint(tmpPoint);

     tmpPoint.setRA(1);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(300.5);
     tmpTrack.addPoint(tmpPoint);
     allTracks.push_back(tmpTrack);

     std::vector<std::pair <unsigned int, unsigned int> > pairs =
	  fieldProximity(allTracks, queryFields, 0);
     
     BOOST_CHECK(pairs.size() == 1);
     BOOST_CHECK(containsPair(1,42,pairs));
     
}

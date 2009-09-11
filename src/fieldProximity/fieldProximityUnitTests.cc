#define BOOST_TEST_MODULE fieldProximityTests

#include <boost/test/included/unit_test.hpp>
#include <boost/current_function.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>


#include "../Exceptions.h"
#include "Field.h"
#include "fieldProximity.h"
#include "TrackForFieldProximity.h"


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
     tmpField.setFieldID("01");
     tmpField.setEpochMJD(300);
     tmpField.setRA(50);
     tmpField.setDec(50);
     tmpField.setRadius(10);
     queryFields.push_back(tmpField);
     
     FieldProximityTrack tmpTrack;
     FieldProximityPoint tmpPoint;
     tmpPoint.setRA(30);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(299);
     tmpTrack.addPoint(tmpPoint);

     tmpPoint.setRA(70);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(301);
     tmpTrack.addPoint(tmpPoint);

     allTracks.push_back(tmpTrack);

     std::vector<std::pair <unsigned int, unsigned int> > pairs =
	  fieldProximity(allTracks, queryFields, 0);
     
     BOOST_CHECK(pairs.size() == 1);
     BOOST_CHECK(containsPair(0,0,pairs));
     
}




BOOST_AUTO_TEST_CASE ( fieldProximity2 ) 
{
     //simple test - one field, two tracks.
     std::vector<Field> queryFields;
     std::vector<FieldProximityTrack> allTracks;

     Field tmpField;
     tmpField.setFieldID("01");
     tmpField.setEpochMJD(300);
     tmpField.setRA(50);
     tmpField.setDec(50);
     tmpField.setRadius(10);
     queryFields.push_back(tmpField);
     
     FieldProximityTrack tmpTrack;
     FieldProximityPoint tmpPoint;
     tmpPoint.setRA(30);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(299);
     tmpTrack.addPoint(tmpPoint);

     tmpPoint.setRA(70);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(301);
     tmpTrack.addPoint(tmpPoint);
     allTracks.push_back(tmpTrack);
     
     //note that this one is too far north
     FieldProximityTrack tmpTrack2;
     tmpPoint.setRA(30);
     tmpPoint.setDec(70);
     tmpPoint.setEpochMJD(299);
     tmpTrack2.addPoint(tmpPoint);

     tmpPoint.setRA(70);
     tmpPoint.setDec(80);
     tmpPoint.setEpochMJD(301);
     tmpTrack2.addPoint(tmpPoint);
     allTracks.push_back(tmpTrack2);     

     std::vector<std::pair <unsigned int, unsigned int> > pairs =
	  fieldProximity(allTracks, queryFields, 0);
     
     BOOST_CHECK(pairs.size() == 1);
     BOOST_CHECK(containsPair(0,0,pairs));
     
}





BOOST_AUTO_TEST_CASE ( fieldProximity3 ) 
{
     //simple test - two fields, 4 tracks.
     std::vector<Field> queryFields;
     std::vector<FieldProximityTrack> allTracks;

     Field tmpField;
     tmpField.setFieldID("01");
     tmpField.setEpochMJD(300);
     tmpField.setRA(50);
     tmpField.setDec(50);
     tmpField.setRadius(10);
     queryFields.push_back(tmpField);

     tmpField.setFieldID("02");
     tmpField.setEpochMJD(300);
     tmpField.setRA(10);
     tmpField.setDec(-30);
     tmpField.setRadius(10);
     queryFields.push_back(tmpField);

     
     FieldProximityTrack tmpTrack;
     FieldProximityPoint tmpPoint;
     tmpPoint.setRA(30);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(299);
     tmpTrack.addPoint(tmpPoint);

     tmpPoint.setRA(70);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(301);
     tmpTrack.addPoint(tmpPoint);
     allTracks.push_back(tmpTrack);
     
     //note that this one is too far north
     FieldProximityTrack tmpTrack2;
     tmpPoint.setRA(30);
     tmpPoint.setDec(70);
     tmpPoint.setEpochMJD(299);
     tmpTrack2.addPoint(tmpPoint);

     tmpPoint.setRA(70);
     tmpPoint.setDec(80);
     tmpPoint.setEpochMJD(301);
     tmpTrack2.addPoint(tmpPoint);
     allTracks.push_back(tmpTrack2);     

     FieldProximityTrack tmpTrack3;
     tmpPoint.setRA(10);
     tmpPoint.setDec(-30);
     tmpPoint.setEpochMJD(299);
     tmpTrack3.addPoint(tmpPoint);

     tmpPoint.setRA(10.000001);
     tmpPoint.setDec(-30.000001);
     tmpPoint.setEpochMJD(301);
     tmpTrack3.addPoint(tmpPoint);
     allTracks.push_back(tmpTrack3);
     
     //note that this one is too far east
     FieldProximityTrack tmpTrack4;
     tmpPoint.setRA(30);
     tmpPoint.setDec(-10);
     tmpPoint.setEpochMJD(299);
     tmpTrack4.addPoint(tmpPoint);

     tmpPoint.setRA(30.5);
     tmpPoint.setDec(-10);
     tmpPoint.setEpochMJD(301);
     tmpTrack4.addPoint(tmpPoint);
     allTracks.push_back(tmpTrack4);     


     std::vector<std::pair <unsigned int, unsigned int> > pairs =
	  fieldProximity(allTracks, queryFields, 0);
     
     BOOST_CHECK(pairs.size() == 2);
     BOOST_CHECK(containsPair(0,0,pairs));
     BOOST_CHECK(containsPair(1,2,pairs));
}





BOOST_AUTO_TEST_CASE ( fieldProximity4 ) 
{
     //try to cause trouble - give no tracks.
     std::vector<Field> queryFields;
     std::vector<FieldProximityTrack> allTracks;

     Field tmpField;
     tmpField.setFieldID("01");
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
     tmpPoint.setEpochMJD(299);
     tmpTrack.addPoint(tmpPoint);

     tmpPoint.setRA(70);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(301);
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
     tmpField.setFieldID("01");
     tmpField.setEpochMJD(300);
     tmpField.setRA(10);
     tmpField.setDec(50);
     tmpField.setRadius(10);
     queryFields.push_back(tmpField);
     
     FieldProximityTrack tmpTrack;
     FieldProximityPoint tmpPoint;
     tmpPoint.setRA(357);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(298);
     tmpTrack.addPoint(tmpPoint);

     tmpPoint.setRA(359);
     tmpPoint.setDec(50);
     tmpPoint.setEpochMJD(299);
     tmpTrack.addPoint(tmpPoint);

     // the next day it will reach position 1, 50
     allTracks.push_back(tmpTrack);

     std::vector<std::pair <unsigned int, unsigned int> > pairs =
	  fieldProximity(allTracks, queryFields, 0);
     
     BOOST_CHECK(pairs.size() == 1);
     BOOST_CHECK(containsPair(0,0,pairs));
     
}

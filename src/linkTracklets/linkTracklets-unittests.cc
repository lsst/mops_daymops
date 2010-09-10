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

// for rand()
#include <cstdlib> 
// for printing timing info
#include <time.h>




#include "lsst/mops/TrackSet.h"
#include "lsst/mops/MopsDetection.h"
#include "lsst/mops/Tracklet.h"
#include "lsst/mops/daymops/linkTracklets/linkTracklets.h"
#include "lsst/mops/Exceptions.h"


namespace lsst {
    namespace mops {



bool Eq(double a, double b) 
{
    double epsilon = 1e-10;
    return (fabs(a - b) < epsilon);
}


void debugPrintTrackletsAndDets(std::vector<MopsDetection> allDets, std::vector<Tracklet> allTracklets) 
{

    for (unsigned int i = 0; i < allDets.size(); i++) {
        MopsDetection* curDet = &allDets.at(i);
        std::cout << curDet->getID() << "\t" << curDet->getRA() << "\t" << curDet->getDec() << '\n';
    }
    std::cout << "all tracklets:\n";

    for (unsigned int i = 0; i < allTracklets.size(); i++) {
        Tracklet* curTracklet = &allTracklets.at(i);
        std::set<unsigned int>::const_iterator dIter;
        for (dIter = curTracklet->indices.begin(); dIter != curTracklet->indices.end(); dIter++) {
            std::cout << *dIter << " ";
        }
        std::cout << '\n';
    }

}

void debugPrintTrackSet(const TrackSet &tracks, const std::vector<MopsDetection> &allDets) 
{
    std::set<Track>::const_iterator trackIter;
    unsigned int trackCount = 0;
    for (trackIter = tracks.componentTracks.begin();
         trackIter != tracks.componentTracks.end();
         trackIter++) {
        std::cout << " track " << trackCount << ":\n";
        std::set<unsigned int>::const_iterator detIdIter;
        for (detIdIter = trackIter->componentDetectionIndices.begin();
             detIdIter != trackIter->componentDetectionIndices.end();
             detIdIter++) {
            std::cout << '\t' << *detIdIter << ": " << allDets.at(*detIdIter).getID() << " "
                      << allDets.at(*detIdIter).getEpochMJD() << 
                " " << allDets.at(*detIdIter).getRA() << " "<< allDets.at(*detIdIter).getDec() << std::endl;
        }
        trackCount++;
    }

}











/*
 * generateTrack: create a new ground-truth track with specified
 * sky-plane behavior and specified obstimes.
 *
 * for each inner vector of trackletObsTimes, we will create a tracklet
 * and add it to allTracklets. It will be given detections at times specified
 * by the doubles contained in the vector.
 * 
 * we will ASSUME that the track is expected to have RA0, Dec0 at time
 * trackletObsTimes[0][0].  If trackletObsTimes is in a non-sorted order, expect
 * weird behavior.
 *
 * each detection will be given a new, unique ID > lastDetId and lastDetId will
 * be MODIFIED to be the last, greated detection ID created. allMopsDetections will
 * be MODIFIED and the new detection will be added to it.
 *
 * similarly, allTracklets will be MODIFIED with new tracklets. Each tracklet
 * will have a new, unique ID > lastTrackletID, and lastTrackletId will be
 * MODIFIED to the last, largest ID we assign.
 *
 * Finally, return a Track containing the IDs of all detections and tracklets
 * created.
 * 
 */
Track generateTrack(double ra0, double dec0, double raV, double decV,
                    double raAcc, double decAcc,
                    std::vector<std::vector <double> > trackletObsTimes,
                    std::vector<MopsDetection> &allDetections,
                    std::vector<Tracklet> &allTracklets, 
                    unsigned int & lastDetId,
                    unsigned int & lastTrackletId) {

    if (trackletObsTimes.size() == 0) {
        throw LSST_EXCEPT(BadParameterException, 
                          std::string(__FUNCTION__)+
                          std::string(": cannot build a track with 0 obs times!"));
    }
    if (trackletObsTimes.at(0).size() == 0) {
        throw LSST_EXCEPT(BadParameterException, 
                          std::string(__FUNCTION__)+
                          std::string(": cannot build a tracklet with 0 obs times!"));
    }
    
    double time0 = trackletObsTimes.at(0).at(0);
    std::vector<std::vector <double > >::const_iterator trackletIter;
    std::vector<double>::const_iterator obsTime;

    Track newTrack;
    for(trackletIter = trackletObsTimes.begin();
        trackletIter != trackletObsTimes.end();
        trackletIter++) {
        if (trackletIter->size() == 0) {
            throw LSST_EXCEPT(BadParameterException, 
                              std::string(__FUNCTION__)+
                              std::string(": cannot build tracklet with 0 obs times!"));
        }
        lastTrackletId++;
        Tracklet newTracklet;
        for (obsTime = trackletIter->begin();
             obsTime != trackletIter->end();
             obsTime++) {
            double resultRa = ra0;
            double resultDec = dec0;
            double tempRaV = raV;
            double tempDecV = decV;
            // calculate ra, dec at this obs time.
            double deltaTime = *obsTime - time0;
            modifyWithAcceleration(resultRa,  tempRaV,  raAcc,  deltaTime);
            modifyWithAcceleration(resultDec, tempDecV, decAcc, deltaTime);
            // create new det, add it to our total set of dets,
            // add its ID to the cur tracklet, and cur track.
            lastDetId++;
            MopsDetection newDet(lastDetId, *obsTime, resultRa, resultDec);
            allDetections.push_back(newDet);
            newTracklet.indices.insert(lastDetId);
            newTrack.componentDetectionIndices.insert(lastDetId);           
        }
        allTracklets.push_back(newTracklet);
        if (allTracklets.size() -1 != lastTrackletId) {
            throw LSST_EXCEPT(BadParameterException,
                              std::string(__FUNCTION__)+
                              std::string(": tracklet IDs are assumed to be the index of the tracklet into the tracklet vector."));
        }
        newTrack.componentTrackletIndices.insert(lastTrackletId);
        
    }

    return newTrack;
}





BOOST_AUTO_TEST_CASE( linkTracklets_realworld_1) {

    std::vector<MopsDetection> allDetections;
    std::vector<Tracklet> allTracklets;
    MopsDetection tmp;
    Tracklet tmpTracklet;

    /*
      2 49523.416229 353.624260 0.521654 20.940280 566 10682481 0. 0.
      74408 49523.422041 353.624302 0.521646 20.352173 566 10682481 0. 0.
      97273 49534.303363 352.326802 1.872652 20.845180 566 10752581 0. 0.
      99300 49534.305146 352.326796 1.872652 20.845182 566 10752581 0. 0.
      130073 49543.427348 350.757850 4.674333 21.576199 566 6206408 0. 0.
      132286 49543.440873 350.755445 4.679549 20.568481 566 2079498 0. 0.
      
      the above track *IS* good enough (though incorrect) but it doesn't get found.
      let's see why.
     */

    linkTrackletsConfig myConfig;
    myConfig.velocityErrorThresh = .192;
    myConfig.detectionLocationErrorThresh = .002;

    tmp.fromMITIString("2 49523.416229 353.624260 0.521654 20.940280 566 10682481 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("74408 49523.422041 353.624302 0.521646 20.352173 566 10682481 0. 0.");
    allDetections.push_back(tmp);

    tmpTracklet.indices.insert(0);
    tmpTracklet.indices.insert(1);
    allTracklets.push_back(tmpTracklet);
    tmpTracklet.indices.clear();

    tmp.fromMITIString("97273 49534.303363 352.326802 1.872652 20.845180 566 10752581 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("99300 49534.305146 352.326796 1.872652 20.845182 566 10752581 0. 0.");
    allDetections.push_back(tmp);

    tmpTracklet.indices.insert(2);
    tmpTracklet.indices.insert(3);
    allTracklets.push_back(tmpTracklet);
    tmpTracklet.indices.clear();

    tmp.fromMITIString("130073 49543.427348 350.757850 4.674333 21.576199 566 6206408 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("132286 49543.440873 350.755445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);

    tmpTracklet.indices.insert(4);
    tmpTracklet.indices.insert(5);
    allTracklets.push_back(tmpTracklet);
    tmpTracklet.indices.clear();

    TrackSet * results = linkTracklets(allDetections, allTracklets, myConfig);
    BOOST_CHECK(results->size() == 1);
   
}




BOOST_AUTO_TEST_CASE( track_1) {
    Track t1;
    Track t2;
    t1.componentDetectionIndices.insert(1);
    t2.componentDetectionIndices.insert(1);
    t1.componentDetectionIndices.insert(2);
    t2.componentDetectionIndices.insert(2);

    t1.componentTrackletIndices.insert(1);
    t2.componentTrackletIndices.insert(1);

    BOOST_CHECK( t1 == t2 );
}


BOOST_AUTO_TEST_CASE( track_2) {
    Track t1;
    Track t2;
    t1.componentDetectionIndices.insert(1);
    t2.componentDetectionIndices.insert(1);
    t1.componentDetectionIndices.insert(2);
    t2.componentDetectionIndices.insert(3);

    t1.componentTrackletIndices.insert(1);
    t2.componentTrackletIndices.insert(2);

    BOOST_CHECK( t1 != t2 );
}


BOOST_AUTO_TEST_CASE( trackSet_1) {
    Track t1;
    Track t2;
    t1.componentDetectionIndices.insert(1);
    t2.componentDetectionIndices.insert(1);
    t1.componentDetectionIndices.insert(2);
    t2.componentDetectionIndices.insert(2);

    t1.componentTrackletIndices.insert(1);
    t2.componentTrackletIndices.insert(1);

    TrackSet ts1;
    TrackSet ts2;
    ts1.insert(t1);
    ts2.insert(t2);

    BOOST_CHECK( ts1 == ts2);
}




BOOST_AUTO_TEST_CASE( trackSet_2) {
    Track t1;
    Track t2;
    t1.componentDetectionIndices.insert(1);
    t2.componentDetectionIndices.insert(1);
    t1.componentDetectionIndices.insert(2);
    t2.componentDetectionIndices.insert(3);

    t1.componentTrackletIndices.insert(1);
    t2.componentTrackletIndices.insert(2);

    TrackSet ts1;
    TrackSet ts2;
    ts1.insert(t1);
    ts2.insert(t2);

    BOOST_CHECK( ts1 != ts2);
}


BOOST_AUTO_TEST_CASE (trackSet_3) {

    Track t1;
    Track t2;
    Track t11;
    Track t22;

    Track t3;

    t1.componentDetectionIndices.insert(1);
    t1.componentDetectionIndices.insert(2);

    t1.componentTrackletIndices.insert(1);


    t2.componentDetectionIndices.insert(3);
    t2.componentDetectionIndices.insert(4);

    t2.componentTrackletIndices.insert(2);


    t11.componentDetectionIndices.insert(10);
    t11.componentDetectionIndices.insert(20);

    t11.componentTrackletIndices.insert(10);


    t22.componentDetectionIndices.insert(30);
    t22.componentDetectionIndices.insert(40);

    t22.componentTrackletIndices.insert(20);


    t3.componentDetectionIndices.insert(4);
    t3.componentDetectionIndices.insert(5);

    t3.componentTrackletIndices.insert(3);

    TrackSet ts1;
    TrackSet ts2;

    ts1.insert(t1);
    ts1.insert(t2);
    BOOST_CHECK(ts1.size() == 2);
    ts2.insert(t1);
    ts2.insert(t2);
    BOOST_CHECK(ts2.size() == 2);
    
    BOOST_CHECK(ts1 == ts2);
    
    ts2.insert(t3);
    BOOST_CHECK(ts2.size() == 3);

    BOOST_CHECK(ts1 != ts2);
    BOOST_CHECK(ts1.isSubsetOf(ts2));

}





BOOST_AUTO_TEST_CASE( linkTracklets_whitebox_getBestFitVelocityAndAcceleration_test0) 
{
    std::vector<double> positions;
    std::vector<double> times;
    positions.push_back(0);
    positions.push_back(1.5);
    positions.push_back(4);
    times.push_back(0);
    times.push_back(1);
    times.push_back(2);
    double velocity, acceleration, position0;
    getBestFitVelocityAndAcceleration(positions, times, velocity, acceleration, position0);
    //std::cout << "position = " << position0 << " + " << velocity << "*t + .5*" << acceleration << "*t^2" << std::endl;
    BOOST_CHECK(Eq(position0,    0));
    BOOST_CHECK(Eq(velocity,     1));
    BOOST_CHECK(Eq(acceleration, 1));    
}





BOOST_AUTO_TEST_CASE( linkTracklets_whitebox_getBestFitVelocityAndAcceleration_test1) 
{
    std::vector<double> positions;
    std::vector<double> times;
    positions.push_back(1.5);
    positions.push_back(4);
    positions.push_back(7.5);
    times.push_back(1);
    times.push_back(2);
    times.push_back(3);
    double velocity, acceleration, position0;
    getBestFitVelocityAndAcceleration(positions, times, velocity, acceleration, position0);
    //std::cout << "position = " << position0 << " + " << velocity << "*t + .5*" << acceleration << "*t^2" << std::endl;
    BOOST_CHECK(Eq(position0,    0));
    BOOST_CHECK(Eq(velocity,     1));
    BOOST_CHECK(Eq(acceleration, 1));    
}






// helper function for creating sets of detections
void addDetectionAt(double MJD, double RA, double dec,  std::vector<MopsDetection> &detVec)
{
    MopsDetection tmpDet(detVec.size(), MJD, RA, dec);
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
  std::vector<MopsDetection> myDets;
  std::vector<Tracklet> pairs;
  linkTrackletsConfig myConfig;
  TrackSet * results = linkTracklets(myDets, pairs, myConfig);
  BOOST_CHECK(pairs.size() == 0);
}







BOOST_AUTO_TEST_CASE( linkTracklets_easy_1 )
{
  std::vector<MopsDetection> myDets;
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

  TrackSet * results = linkTracklets(myDets, pairs, myConfig);

  BOOST_CHECK(results->size() == 1);
}







BOOST_AUTO_TEST_CASE( linkTracklets_easy_2 )
{
    // same as 1, but with more tracks (all clearly separated)

  std::vector<MopsDetection> myDets;
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

  TrackSet * results = linkTracklets(myDets, pairs, myConfig);

  BOOST_CHECK(results->size() == 10);
}






BOOST_AUTO_TEST_CASE( linkTracklets_easy_3 )
{
    // same as 2, but with tracks crossing RA 0 line

  std::vector<MopsDetection> myDets;
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

  TrackSet * results = linkTracklets(myDets, pairs, myConfig);

  BOOST_CHECK(results->size() == 10);
}





BOOST_AUTO_TEST_CASE( linkTracklets_easy_4_1 )
{
    // same as 1, but with track crossing RA 0 line

  std::vector<MopsDetection> myDets;
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

  TrackSet * results = linkTracklets(myDets, pairs, myConfig);

  //std::cout << "results size = " << results->size() << '\n';
  BOOST_CHECK(results->size() == 1);
}



BOOST_AUTO_TEST_CASE( linkTracklets_easy_4 )
{
    // same as 2, but with tracks crossing RA 0 line

  std::vector<MopsDetection> myDets;
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

  TrackSet * results = linkTracklets(myDets, pairs, myConfig);

  //std::cout << "results size = " << results->size() << '\n';
  BOOST_CHECK(results->size() == 10);
}



BOOST_AUTO_TEST_CASE( linkTracklets_easy_5 )
{
    // same as 4, but going the other way!

  std::vector<MopsDetection> myDets;
  std::vector<Tracklet> pairs;
  for (unsigned int i = 0; i < 10; i++) {

      addDetectionAt(5300.0,    0.101,    50.201 + i, myDets);
      addDetectionAt(5300.01,   0.1,      50.2   + i, myDets);

      addDetectionAt(5301.0,    0.001,    50.101 + i, myDets);
      addDetectionAt(5301.01,   0.,       50.1   + i, myDets);

      addDetectionAt(5302.0,  359.901,    50.001 + i, myDets);
      addDetectionAt(5302.01, 359.9,      50.    + i, myDets);

      addPair(0 + 6*i,   1 + 6*i,   pairs);
      addPair(2 + 6*i,   3 + 6*i,   pairs);
      addPair(4 + 6*i,   5 + 6*i,   pairs);
      
  }

  
  linkTrackletsConfig myConfig;

  TrackSet * results = linkTracklets(myDets, pairs, myConfig);

  //std::cout << "results size = " << results->size() << '\n';
  // for (unsigned int i = 0; i < results->size(); i++) {
      
  //     std::set<unsigned int>::const_iterator dIter;
  //     const std::set<unsigned int> * cur = &(results->at(i).componentDetectionIndices);

  //     for (dIter = cur->begin(); dIter != cur->end(); dIter++) {
  //         std::cout << " " << *dIter;
  //     }
  //     std::cout << "\n";
  //     cur = & ( results->at(i).componentTrackletIndices);
  //     for (dIter = cur->begin(); dIter != cur->end(); dIter++) {
  //         std::cout << " " << *dIter;
  //     }
  //     std::cout << "\n";

  // }

  BOOST_CHECK(results->size() == 10);
}



// Track generateTrack(double ra0, double dec0, double raV, double decV,
//                     double raAcc, double decAcc,
//                     std::vector<std::vector <double> > trackletObsTimes,
//                     std::vector<MopsDetection> &allMopsDetections,
//                     std::vector<Tracklet> &allTracklets, 
//                     unsigned int & lastDetId,
//                     unsigned int & lastTrackletId) {

BOOST_AUTO_TEST_CASE( linkTracklets_1 )
{
    TrackSet expectedTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;

    std::vector<std::vector<double> > imgTimes(3);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.03);

    imgTimes.at(1).push_back(5301);
    imgTimes.at(1).push_back(5301.03);

    imgTimes.at(2).push_back(5302);
    imgTimes.at(2).push_back(5302.03);

    expectedTracks.insert(generateTrack(20., 20., 
                                        .25, .01, 
                                        .0002, .002,
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));



    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);

    BOOST_CHECK(*foundTracks == expectedTracks);

}


BOOST_AUTO_TEST_CASE( linkTracklets_2 )
{
    // lots of support nodes!
    TrackSet expectedTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;

    std::vector<std::vector<double> > imgTimes(7);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.03);

    imgTimes.at(1).push_back(5301);
    imgTimes.at(1).push_back(5301.03);

    imgTimes.at(2).push_back(5302);
    imgTimes.at(2).push_back(5302.03);

    imgTimes.at(3).push_back(5303);
    imgTimes.at(3).push_back(5303.03);

    imgTimes.at(4).push_back(5304);
    imgTimes.at(4).push_back(5304.03);

    imgTimes.at(5).push_back(5305);
    imgTimes.at(5).push_back(5305.03);

    imgTimes.at(6).push_back(5306);
    imgTimes.at(6).push_back(5306.03);

    expectedTracks.insert(generateTrack(20., 20., 
                                        .25, .01, 
                                        .0002, .002,
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));



    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);


    BOOST_CHECK(expectedTracks.isSubsetOf(*foundTracks));

}



BOOST_AUTO_TEST_CASE( linkTracklets_3 )
{
    // lots of support nodes, and a psuedo-deep stack
    TrackSet expectedTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;

    std::vector<std::vector<double> > imgTimes(7);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.03);

    imgTimes.at(1).push_back(5301);
    imgTimes.at(1).push_back(5301.03);

    imgTimes.at(2).push_back(5302);
    imgTimes.at(2).push_back(5302.03);

    imgTimes.at(3).push_back(5303);
    imgTimes.at(3).push_back(5303.005);    
    imgTimes.at(3).push_back(5303.01);
    imgTimes.at(3).push_back(5303.015);
    imgTimes.at(3).push_back(5303.02);
    imgTimes.at(3).push_back(5303.025);
    imgTimes.at(3).push_back(5303.03);
    imgTimes.at(3).push_back(5303.035);
    imgTimes.at(3).push_back(5303.04);
    imgTimes.at(3).push_back(5303.045);
    imgTimes.at(3).push_back(5303.05);
    imgTimes.at(3).push_back(5303.055);

    imgTimes.at(4).push_back(5304);
    imgTimes.at(4).push_back(5304.03);

    imgTimes.at(5).push_back(5305);
    imgTimes.at(5).push_back(5305.03);

    imgTimes.at(6).push_back(5306);
    imgTimes.at(6).push_back(5306.03);

    expectedTracks.insert(generateTrack(20., 20., 
                                        .25, .01, 
                                        .0002, .002,
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));



    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);

    

    BOOST_CHECK(expectedTracks.isSubsetOf(*foundTracks));
    // std::cout << " got tracks: \n";
    // debugPrintTrackSet(foundTracks, allDets);
    // std::cout << " expected tracks: \n";
    // debugPrintTrackSet(expectedTracks, allDets);

}



BOOST_AUTO_TEST_CASE( linkTracklets_4 )
{

    // lots of tracks this time. still simple cadence.
    TrackSet expectedTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;

    std::vector<std::vector<double> > imgTimes(3);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.03);

    imgTimes.at(1).push_back(5301);
    imgTimes.at(1).push_back(5301.03);

    imgTimes.at(2).push_back(5302);
    imgTimes.at(2).push_back(5302.03);

    expectedTracks.insert(generateTrack(20., 20., 
                                        .25, .01, 
                                        .0002, .002,
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));

    expectedTracks.insert(generateTrack(10., 10., 
                                        -.02, .015, 
                                        .00025, .0002,
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));

    expectedTracks.insert(generateTrack(15., 10., 
                                        .001, .00001, 
                                        .00015, .000023,
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));

    expectedTracks.insert(generateTrack(12.5, 12.5, 
                                        -.01, -.001, 
                                        .000001, -.00023,
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));

    expectedTracks.insert(generateTrack(13, 14., 
                                        .001, -.01, 
                                        -.001, -.00023,
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));

    expectedTracks.insert(generateTrack(14.5, 14., 
                                        .01, -.000001, 
                                        -.00015, .00023,
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));

    expectedTracks.insert(generateTrack(16.5, 14., 
                                        .0155, .000001, 
                                        -.00015, -.00023,
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));

    expectedTracks.insert(generateTrack(10.5, 13., 
                                        .000155, .0001, 
                                        .00066, -.00066,
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));

    expectedTracks.insert(generateTrack(17.5, 12.5, 
                                        -.011, -.001, 
                                        .0001112, -.0002388,
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));


    expectedTracks.insert(generateTrack(12, 19.5, 
                                        .001333, -.008888, 
                                        -.0039083, -.001999,
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));

    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);

    BOOST_CHECK(*foundTracks == expectedTracks);

}







BOOST_AUTO_TEST_CASE( linkTracklets_4_pt_5 )
{

    // lots of tracks, following a coherent pattern but randomly perturbed.
    TrackSet expectedTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;

    std::vector<std::vector<double> > imgTimes(3);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.03);

    imgTimes.at(1).push_back(5305);
    imgTimes.at(1).push_back(5305.03);

    imgTimes.at(2).push_back(5312);
    imgTimes.at(2).push_back(5312.03);

    // seed the random number generator with a known value;
    // this way the test will be identical on each run.
    srand(2); // srand(1);

    for (unsigned int i = 0; i < 10; i++) {

        // get 6 floating point numbers between 0 and 1.
        std::vector<double> someRands;         
        for (unsigned int j = 0; j < 6; j++) {
            someRands.push_back( (double) rand() / RAND_MAX );
        }
        //generate a random permutation on this track.
        expectedTracks.insert(generateTrack(20. + someRands[0] * 10., //location 
                                            20. + someRands[1] * 10., //location
                                            (someRands[2] - .1) * 2., //1 in 10 chance of retrograde, maxv 2
                                            (someRands[3] - .5) * .5, // maxv .5 in dec 
                                            (someRands[4]) * .0019, //max acc of .0019, always positive
                                            (someRands[5]) * .0019, //same
                                            imgTimes, 
                                            allDets, allTracklets, 
                                            firstDetId, firstTrackletId));
    }

    struct tm * timeinfo;
    time_t rawtime;
    
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "Generated " << expectedTracks.size() << " tracks. Calling linkTracklets ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "got " << foundTracks->size() << " results, checking if they contain the true tracks ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    BOOST_CHECK(expectedTracks.isSubsetOf(*foundTracks));

    // std::cout << " got tracks: \n";
    // debugPrintTrackSet(foundTracks, allDets);
    // std::cout << " expected tracks: \n";
    // debugPrintTrackSet(expectedTracks, allDets);



}





BOOST_AUTO_TEST_CASE( linkTracklets_5 )
{

    // " << expectedTracks.size() << " tracks, following a coherent pattern but randomly perturbed.
    TrackSet expectedTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;

    std::vector<std::vector<double> > imgTimes(3);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.03);

    imgTimes.at(1).push_back(5305);
    imgTimes.at(1).push_back(5305.03);

    imgTimes.at(2).push_back(5312);
    imgTimes.at(2).push_back(5312.03);

    // seed the random number generator with a known value;
    // this way the test will be identical on each run.
    srand(1);

    for (unsigned int i = 0; i < 100; i++) {

        // get 6 floating point numbers between 0 and 1.
        std::vector<double> someRands;         
        for (unsigned int j = 0; j < 6; j++) {
            someRands.push_back( (double) rand() / RAND_MAX );
        }
        //generate a random permutation on this track.
        expectedTracks.insert(generateTrack(20. + someRands[0] * 10., //location 
                                            20. + someRands[1] * 10., //location
                                            (someRands[2] - .1) * 2., //1 in 10 chance of retrograde, maxv 2
                                            (someRands[3] - .5) * .5, // maxv .5 in dec 
                                            (someRands[4]) * .0019, //max acc of .0019, always positive
                                            (someRands[5]) * .0019, //same
                                            imgTimes, 
                                            allDets, allTracklets, 
                                            firstDetId, firstTrackletId));

    }

    struct tm * timeinfo;
    time_t rawtime;
    
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "Generated " << expectedTracks.size() << " tracks. Calling linkTracklets ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "got " << foundTracks->size() << " results, checking if they contain the true tracks ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    BOOST_CHECK(expectedTracks.isSubsetOf(*foundTracks));

}





BOOST_AUTO_TEST_CASE( linkTracklets_5_high_acc )
{

    // " << expectedTracks.size() << " tracks, following a coherent pattern but randomly perturbed.
    TrackSet expectedTracks;
    TrackSet implausibleTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;
    myConfig.leafSize=16;

    myConfig.maxRAAccel = .02;
    myConfig.maxDecAccel = .02;
    
    std::vector<std::vector<double> > imgTimes(3);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.01);

    imgTimes.at(1).push_back(5315);
    imgTimes.at(1).push_back(5315.01);

    imgTimes.at(2).push_back(5330);
    imgTimes.at(2).push_back(5330.01);


    expectedTracks.insert(generateTrack(20., //location ra 
                                        20., //location dec
                                        .05, // v0 ra
                                        .05, // v0 dec 
                                        .0199, //ra acc
                                        .0199, //dec acc
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));
    
    implausibleTracks.insert(generateTrack(10.,
                                           10.,
                                           .05,
                                           .05,
                                           .021,
                                           .021,
                                           imgTimes,
                                           allDets, allTracklets,
                                           firstDetId, firstTrackletId));

    
    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);

    BOOST_CHECK(expectedTracks.isSubsetOf(*foundTracks));

    BOOST_CHECK(!implausibleTracks.isSubsetOf(*foundTracks));

}





BOOST_AUTO_TEST_CASE( linkTracklets_too_high_acc )
{

    // " << expectedTracks.size() << " tracks, following a coherent pattern but randomly perturbed.
    TrackSet implausibleTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;

    myConfig.maxRAAccel = .02;
    myConfig.maxDecAccel = .02;
    
    std::vector<std::vector<double> > imgTimes(3);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.01);

    imgTimes.at(1).push_back(5315);
    imgTimes.at(1).push_back(5315.01);

    imgTimes.at(2).push_back(5330);
    imgTimes.at(2).push_back(5330.01);


    
    implausibleTracks.insert(generateTrack(10.,
                                           10.,
                                           .05,
                                           .05,
                                           .021,
                                           .021,
                                           imgTimes,
                                           allDets, allTracklets,
                                           firstDetId, firstTrackletId));

    
    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);


    BOOST_CHECK(!implausibleTracks.isSubsetOf(*foundTracks));

}





BOOST_AUTO_TEST_CASE( linkTracklets_5_1_plus_high_acc )
{

    // " << expectedTracks.size() << " tracks, following a coherent pattern but randomly perturbed.
    TrackSet expectedTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;

    std::vector<std::vector<double> > imgTimes(3);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.03);

    imgTimes.at(1).push_back(5305);
    imgTimes.at(1).push_back(5305.03);

    imgTimes.at(2).push_back(5312);
    imgTimes.at(2).push_back(5312.03);

    myConfig.leafSize = 1;

    // seed the random number generator with a known value;
    // this way the test will be identical on each run.
    srand(2);

    for (unsigned int i = 0; i < 200; i++) {

        // get 6 floating point numbers between 0 and 1.
        std::vector<double> someRands;         
        for (unsigned int j = 0; j < 6; j++) {
            someRands.push_back( (double) rand() / RAND_MAX );
        }
        //generate a random permutation on this track.
        expectedTracks.insert(generateTrack(20. + someRands[0] * 10., //location 
                                            20. + someRands[1] * 10., //location
                                            (someRands[2] - .1) * 2., //1 in 10 chance of retrograde, maxv 2
                                            (someRands[3] - .5) * .5, // maxv .5 in dec 
                                            (someRands[4]) * .0019, //max acc of .0019, always positive
                                            (someRands[5]) * .0019, //same
                                            imgTimes, 
                                            allDets, allTracklets, 
                                            firstDetId, firstTrackletId));

    }

    // add a fast mover, too
    expectedTracks.insert(generateTrack(20., //location ra 
                                        20., //location dec
                                        .05, // v0 ra
                                        .05, // v0 dec 
                                        .0199, //ra acc
                                        .0199, //dec acc
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));

    // add a fast mover, too
    expectedTracks.insert(generateTrack(20., //location ra 
                                        20., //location dec
                                        .05, // v0 ra
                                        .05, // v0 dec 
                                        -.0199, //ra acc
                                        -.0199, //dec acc
                                        imgTimes, 
                                        allDets, allTracklets, 
                                        firstDetId, firstTrackletId));

    struct tm * timeinfo;
    time_t rawtime;
    
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "Generated " << expectedTracks.size() << " tracks. Calling linkTracklets ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "got " << foundTracks->size() << " results, checking if they contain the true tracks ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    BOOST_CHECK(expectedTracks.isSubsetOf(*foundTracks));

}






BOOST_AUTO_TEST_CASE( linkTracklets_5_1 )
{

    // " << expectedTracks.size() << " tracks, following a coherent pattern but randomly perturbed.
    TrackSet expectedTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;

    std::vector<std::vector<double> > imgTimes(3);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.03);

    imgTimes.at(1).push_back(5305);
    imgTimes.at(1).push_back(5305.03);

    imgTimes.at(2).push_back(5312);
    imgTimes.at(2).push_back(5312.03);

    // seed the random number generator with a known value;
    // this way the test will be identical on each run.
    srand(2);

    for (unsigned int i = 0; i < 200; i++) {

        // get 6 floating point numbers between 0 and 1.
        std::vector<double> someRands;         
        for (unsigned int j = 0; j < 6; j++) {
            someRands.push_back( (double) rand() / RAND_MAX );
        }
        //generate a random permutation on this track.
        expectedTracks.insert(generateTrack(20. + someRands[0] * 10., //location 
                                            20. + someRands[1] * 10., //location
                                            (someRands[2] - .1) * 2., //1 in 10 chance of retrograde, maxv 2
                                            (someRands[3] - .5) * .5, // maxv .5 in dec 
                                            (someRands[4]) * .0019, //max acc of .0019, always positive
                                            (someRands[5]) * .0019, //same
                                            imgTimes, 
                                            allDets, allTracklets, 
                                            firstDetId, firstTrackletId));

    }

    struct tm * timeinfo;
    time_t rawtime;
    
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "Generated " << expectedTracks.size() << " tracks. Calling linkTracklets ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "got " << foundTracks->size() << " results, checking if they contain the true tracks ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    BOOST_CHECK(expectedTracks.isSubsetOf(*foundTracks));

}




BOOST_AUTO_TEST_CASE( linkTracklets_5_2 )
{

    // " << expectedTracks.size() << " tracks, following a coherent pattern but randomly perturbed.
    TrackSet expectedTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;

    std::vector<std::vector<double> > imgTimes(3);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.03);

    imgTimes.at(1).push_back(5305);
    imgTimes.at(1).push_back(5305.03);

    imgTimes.at(2).push_back(5312);
    imgTimes.at(2).push_back(5312.03);

    // seed the random number generator with a known value;
    // this way the test will be identical on each run.
    srand(3);

    for (unsigned int i = 0; i < 300; i++) {

        // get 6 floating point numbers between 0 and 1.
        std::vector<double> someRands;         
        for (unsigned int j = 0; j < 6; j++) {
            someRands.push_back( (double) rand() / RAND_MAX );
        }
        //generate a random permutation on this track.
        expectedTracks.insert(generateTrack(20. + someRands[0] * 10., //location 
                                            20. + someRands[1] * 10., //location
                                            (someRands[2] - .1) * 2., //1 in 10 chance of retrograde, maxv 2
                                            (someRands[3] - .5) * .5, // maxv .5 in dec 
                                            (someRands[4]) * .0019, //max acc of .0019, always positive
                                            (someRands[5]) * .0019, //same
                                            imgTimes, 
                                            allDets, allTracklets, 
                                            firstDetId, firstTrackletId));

    }

    struct tm * timeinfo;
    time_t rawtime;
    
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "Generated " << expectedTracks.size() << " tracks. Calling linkTracklets ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "got " << foundTracks->size() << " results, checking if they contain the true tracks ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    BOOST_CHECK(expectedTracks.isSubsetOf(*foundTracks));

}





BOOST_AUTO_TEST_CASE( linkTracklets_5_3 )
{

    // " << expectedTracks.size() << " tracks, following a coherent pattern but randomly perturbed.
    TrackSet expectedTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;

    std::vector<std::vector<double> > imgTimes(3);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.03);

    imgTimes.at(1).push_back(5305);
    imgTimes.at(1).push_back(5305.03);

    imgTimes.at(2).push_back(5312);
    imgTimes.at(2).push_back(5312.03);

    // seed the random number generator with a known value;
    // this way the test will be identical on each run.
    srand(4);

    for (unsigned int i = 0; i < 400; i++) {

        // get 6 floating point numbers between 0 and 1.
        std::vector<double> someRands;         
        for (unsigned int j = 0; j < 6; j++) {
            someRands.push_back( (double) rand() / RAND_MAX );
        }
        //generate a random permutation on this track.
        expectedTracks.insert(generateTrack(20. + someRands[0] * 10., //location 
                                            20. + someRands[1] * 10., //location
                                            (someRands[2] - .1) * 2., //1 in 10 chance of retrograde, maxv 2
                                            (someRands[3] - .5) * .5, // maxv .5 in dec 
                                            (someRands[4]) * .0019, //max acc of .0019, always positive
                                            (someRands[5]) * .0019, //same
                                            imgTimes, 
                                            allDets, allTracklets, 
                                            firstDetId, firstTrackletId));

    }

    struct tm * timeinfo;
    time_t rawtime;
    
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "Generated " << expectedTracks.size() << " tracks. Calling linkTracklets ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "got " << foundTracks->size() << " results, checking if they contain the true tracks ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    BOOST_CHECK(expectedTracks.isSubsetOf(*foundTracks));

}

BOOST_AUTO_TEST_CASE( linkTracklets_5_4 )
{

    // " << expectedTracks.size() << " tracks, following a coherent pattern but randomly perturbed.
    TrackSet expectedTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;

    std::vector<std::vector<double> > imgTimes(3);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.03);

    imgTimes.at(1).push_back(5305);
    imgTimes.at(1).push_back(5305.03);

    imgTimes.at(2).push_back(5312);
    imgTimes.at(2).push_back(5312.03);

    // seed the random number generator with a known value;
    // this way the test will be identical on each run.
    srand(5);

    for (unsigned int i = 0; i < 500; i++) {

        // get 6 floating point numbers between 0 and 1.
        std::vector<double> someRands;         
        for (unsigned int j = 0; j < 6; j++) {
            someRands.push_back( (double) rand() / RAND_MAX );
        }
        //generate a random permutation on this track.
        expectedTracks.insert(generateTrack(20. + someRands[0] * 10., //location 
                                            20. + someRands[1] * 10., //location
                                            (someRands[2] - .1) * 2., //1 in 10 chance of retrograde, maxv 2
                                            (someRands[3] - .5) * .5, // maxv .5 in dec 
                                            (someRands[4]) * .0019, //max acc of .0019, always positive
                                            (someRands[5]) * .0019, //same
                                            imgTimes, 
                                            allDets, allTracklets, 
                                            firstDetId, firstTrackletId));

    }

    struct tm * timeinfo;
    time_t rawtime;
    
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "Generated " << expectedTracks.size() << " tracks. Calling linkTracklets ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "got " << foundTracks->size() << " results, checking if they contain the true tracks ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    BOOST_CHECK(expectedTracks.isSubsetOf(*foundTracks));

}






BOOST_AUTO_TEST_CASE( linkTracklets_5_5 )
{

    // " << expectedTracks.size() << " tracks, following a coherent pattern but randomly perturbed.
    TrackSet expectedTracks;
    std::vector<MopsDetection> allDets;
    std::vector<Tracklet> allTracklets;
    unsigned int firstDetId = -1;
    unsigned int firstTrackletId = -1;

  
    linkTrackletsConfig myConfig;

    std::vector<std::vector<double> > imgTimes(3);

    imgTimes.at(0).push_back(5300);
    imgTimes.at(0).push_back(5300.03);

    imgTimes.at(1).push_back(5305);
    imgTimes.at(1).push_back(5305.03);

    imgTimes.at(2).push_back(5312);
    imgTimes.at(2).push_back(5312.03);

    // seed the random number generator with a known value;
    // this way the test will be identical on each run.
    srand(6);

    for (unsigned int i = 0; i < 600; i++) {

        // get 6 floating point numbers between 0 and 1.
        std::vector<double> someRands;         
        for (unsigned int j = 0; j < 6; j++) {
            someRands.push_back( (double) rand() / RAND_MAX );
        }
        //generate a random permutation on this track.
        expectedTracks.insert(generateTrack(20. + someRands[0] * 10., //location 
                                            20. + someRands[1] * 10., //location
                                            (someRands[2] - .1) * 2., //1 in 10 chance of retrograde, maxv 2
                                            (someRands[3] - .5) * .5, // maxv .5 in dec 
                                            (someRands[4]) * .0019, //max acc of .0019, always positive
                                            (someRands[5]) * .0019, //same
                                            imgTimes, 
                                            allDets, allTracklets, 
                                            firstDetId, firstTrackletId));

    }

    struct tm * timeinfo;
    time_t rawtime;
    
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "Generated " << expectedTracks.size() << " tracks. Calling linkTracklets ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    TrackSet * foundTracks = linkTracklets(allDets, allTracklets, myConfig);

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );                    
    
    std::cout << "got " << foundTracks->size() << " results, checking if they contain the true tracks ";
    std::cout << " current wall-clock time is " << asctime (timeinfo) << std::endl;

    BOOST_CHECK(expectedTracks.isSubsetOf(*foundTracks));

}



// TBD: check that tracks with too-high acceleration are correctly rejected, etc.







}} // close lsst::mops

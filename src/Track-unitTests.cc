// -*- LSST-C++ -*-
/* jonathan myers */
#define BOOST_TEST_MODULE Track

#include <boost/test/included/unit_test.hpp>
#include <boost/current_function.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>




#include "lsst/mops/Track.h"
#include "lsst/mops/Tracklet.h"
#include "lsst/mops/TrackSet.h"
#include "lsst/mops/MopsDetection.h"
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
        std::set<unsigned int> componentDetectionIndices = trackIter->getComponentDetectionIndices();
        for (detIdIter = componentDetectionIndices.begin();
             detIdIter != componentDetectionIndices.end();
             detIdIter++) {
            std::cout << '\t' << *detIdIter << ": " << allDets.at(*detIdIter).getID() << " "
                      << allDets.at(*detIdIter).getEpochMJD() << 
                " " << allDets.at(*detIdIter).getRA() << " "<< allDets.at(*detIdIter).getDec() << std::endl;
        }
        trackCount++;
    }

}



BOOST_AUTO_TEST_CASE( track_1) {
    Track t1;
    Track t2;
    std::vector<MopsDetection> allDetections;
    MopsDetection tmp;
    tmp.fromMITIString("130073 49543.427348 350.757850 4.674333 21.576199 566 6206408 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130074 49543.440873 350.755445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130075 49543.440873 350.755445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    t1.addDetection(1, allDetections);
    t2.addDetection(1, allDetections);
    t1.addDetection(2, allDetections);
    t2.addDetection(2, allDetections);


    BOOST_CHECK( t1 == t2 );
}


BOOST_AUTO_TEST_CASE( track_2) {
    Track t1;
    Track t2;
    std::vector<MopsDetection> allDetections;
    MopsDetection tmp;
    tmp.fromMITIString("130073 49543.427348 350.757850 4.674333 21.576199 566 6206408 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130074 49543.440873 350.755445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130075 49543.450873 350.855445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130076 49543.450873 350.855445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);

    t1.addDetection(1, allDetections);
    t2.addDetection(1, allDetections);
    t1.addDetection(2, allDetections);
    t2.addDetection(3, allDetections);

    t1.componentTrackletIndices.insert(1);
    t2.componentTrackletIndices.insert(2);
    BOOST_CHECK( t1 != t2 );
}


BOOST_AUTO_TEST_CASE( trackSet_1) {
    Track t1;
    Track t2;
    std::vector<MopsDetection> allDetections;
    MopsDetection tmp;
    tmp.fromMITIString("130073 49543.427348 350.757850 4.674333 21.576199 566 6206408 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130074 49543.440873 350.755445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130076 49543.450873 350.855445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    t1.addDetection(1, allDetections);
    t2.addDetection(1, allDetections);
    t1.addDetection(2, allDetections);
    t2.addDetection(2, allDetections);

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
    std::vector<MopsDetection> allDetections;
    MopsDetection tmp;
    tmp.fromMITIString("130073 49543.427348 350.757850 4.674333 21.576199 566 6206408 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130074 49543.440873 350.755445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130075 49543.450873 350.855445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130076 49543.450873 350.855445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);

    t1.addDetection(1, allDetections);
    t2.addDetection(1, allDetections);
    t1.addDetection(2, allDetections);
    t2.addDetection(3, allDetections);

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

    MopsDetection tmp;
    std::vector<MopsDetection> allDetections;
    tmp.fromMITIString("130073 49543.427348 350.757850 4.674333 21.576199 566 6206408 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130074 49543.440873 350.755445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130075 49543.450873 350.855445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130076 49544.440873 351.855445 2.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130077 49543.427348 350.757850 4.674333 21.576199 566 6206408 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130078 49543.440873 350.755445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130079 49543.450873 350.855445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130080 49544.440873 351.855445 2.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130081 49543.450873 350.855445 4.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130082 49544.440873 351.855445 2.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130083 49544.440873 351.855445 2.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);
    tmp.fromMITIString("130084 49544.440873 351.855445 2.679549 20.568481 566 2079498 0. 0.");
    allDetections.push_back(tmp);

    t1.addDetection(1, allDetections);
    t1.addDetection(2, allDetections);

    t1.componentTrackletIndices.insert(1);


    t2.addDetection(3, allDetections);
    t2.addDetection(4, allDetections);

    t2.componentTrackletIndices.insert(2);


    t11.addDetection(5, allDetections);
    t11.addDetection(6, allDetections);

    t11.componentTrackletIndices.insert(10);


    t22.addDetection(7, allDetections);
    t22.addDetection(8, allDetections);

    t22.componentTrackletIndices.insert(20);


    t3.addDetection(9, allDetections);
    t3.addDetection(10, allDetections);

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



// helper function for creating sets of detections
void addDetectionAt(double MJD, double RA, double dec,  std::vector<MopsDetection> &detVec)
{
    MopsDetection tmpDet(detVec.size(), MJD, RA, dec);
    detVec.push_back(tmpDet);
}



template <typename T>
bool setsEqual(const std::set<T> s1, const std::set<T> s2)
{
    if (s1.size() != s2.size()) {
        return false;
    }

    // take advantage of the fact that sets are ordered
    typename std::set<T>::const_iterator sIter1;
    typename std::set<T>::const_iterator sIter2;
    sIter1 = s1.begin();
    sIter2 = s2.begin();

    while ((sIter1 != s1.end()) && (sIter2 != s2.end())) {
        if (*sIter1  != *sIter2) {
            return false;
        }
        sIter1++;
        sIter2++;
    }
    return true;
}



BOOST_AUTO_TEST_CASE( track_quadraticFitting_addTracklet ) 
{
    std::vector<MopsDetection> allDets;
    addDetectionAt(0, 0, 0,     allDets);
    addDetectionAt(1, 1.5, 1.5, allDets);
    addDetectionAt(2, 4, 4,     allDets);
    addDetectionAt(3, 4.5, 4.5,     allDets);

    Tracklet t1, t2;
    allDets.at(0).setID(42);
    allDets.at(1).setID(43);
    allDets.at(2).setID(1337);
    allDets.at(3).setID(1338);

    t1.indices.insert(0);
    t1.indices.insert(1);

    t2.indices.insert(2);
    t2.indices.insert(3);

    Track testTrack;

    testTrack.addTracklet(0, t1, allDets);
    testTrack.addTracklet(1, t2, allDets);
    
    std::set<unsigned int> expectedDetIndices;
    expectedDetIndices.insert(0);
    expectedDetIndices.insert(1);
    expectedDetIndices.insert(2);
    expectedDetIndices.insert(3);
    std::set<unsigned int> expectedDiaIds;
    expectedDiaIds.insert(42);
    expectedDiaIds.insert(43);
    expectedDiaIds.insert(1337);
    expectedDiaIds.insert(1338);
    std::set<unsigned int> expectedTrackletIndices;
    expectedTrackletIndices.insert(0);
    expectedTrackletIndices.insert(1);

    BOOST_CHECK(setsEqual(testTrack.componentTrackletIndices, 
                          expectedTrackletIndices));
    
    BOOST_CHECK( setsEqual (testTrack.getComponentDetectionIndices(), 
                            expectedDetIndices));
    
    BOOST_CHECK( setsEqual( testTrack.getComponentDetectionDiaIds(), 
                            expectedDiaIds));
}





// jmyers - these predate the use of topocentric correction in the
// detections, and so they no longer work.  Consider retooling them
// someday!

// BOOST_AUTO_TEST_CASE( track_quadraticFitting_test0) 
// {
//     std::vector<MopsDetection> allDets;
//     addDetectionAt(0, 0, 0,     allDets);
//     addDetectionAt(1, 1.5, 1.5, allDets);
//     addDetectionAt(2, 4, 4,     allDets);
//     Track testTrack;
//     testTrack.addDetection(0, allDets);
//     testTrack.addDetection(1, allDets);
//     testTrack.addDetection(2, allDets);
//     double epoch, ra0, raV, raAcc, dec0, decV, decAcc;
//     testTrack.calculateBestFitQuadratic(allDets);
//     testTrack.getBestFitQuadratic(epoch, ra0, raV, raAcc, dec0, decV, decAcc);
//     BOOST_CHECK(Eq(epoch,   0.));
//     BOOST_CHECK(Eq(ra0,     0.));
//     BOOST_CHECK(Eq(raV,     1.));
//     BOOST_CHECK(Eq(raAcc,   1.));    
//     BOOST_CHECK(Eq(dec0,    0.));
//     BOOST_CHECK(Eq(decV,    1.));
//     BOOST_CHECK(Eq(decAcc,  1.));    
// }





// BOOST_AUTO_TEST_CASE( track_quadraticFit_test1) 
// {
//     std::vector<MopsDetection> allDets;
//     addDetectionAt(1, 1, 1,    allDets);
//     addDetectionAt(2, 2.5, 2.5,allDets);
//     addDetectionAt(3, 5,   5,  allDets);
//     Track testTrack;
//     testTrack.addDetection(0, allDets);
//     testTrack.addDetection(1, allDets);
//     testTrack.addDetection(2, allDets);
//     double epoch, ra0, raV, raAcc, dec0, decV, decAcc;
//     testTrack.calculateBestFitQuadratic(allDets);
//     testTrack.getBestFitQuadratic(epoch, ra0, raV, raAcc, dec0, decV, decAcc);

//     BOOST_CHECK(Eq(epoch,        1));
//     BOOST_CHECK(Eq(ra0,          1));
//     BOOST_CHECK(Eq(raV,          1));
//     BOOST_CHECK(Eq(raAcc,        1));    
//     BOOST_CHECK(Eq(dec0,          1));
//     BOOST_CHECK(Eq(decV,          1));
//     BOOST_CHECK(Eq(decAcc,        1));    
//     double predRa, predDec;
//     testTrack.predictLocationAtTime(1,predRa, predDec);
//     BOOST_CHECK(Eq(predRa,  1));
//     BOOST_CHECK(Eq(predDec, 1));
//     testTrack.predictLocationAtTime(2,predRa, predDec);
//     BOOST_CHECK(Eq(predRa,  2.5));
//     BOOST_CHECK(Eq(predDec, 2.5));
// }



double projectLoc(double time, double p0, double v, double acc)
{
     return p0 + v*time + .5 * time* time * acc;
}

#define NUM_CALLS 100000

BOOST_AUTO_TEST_CASE( big_quadFit_test0) 
{
    // basically copied from r. 16930 daymops/trunk/tests/quadraticFitting/benchmark.cc
    
    std::vector<double> times;
    times.push_back(0.);
    times.push_back(0.05);
    times.push_back(1.);
    times.push_back(1.05);
    times.push_back(5.);
    times.push_back(5.05);
    
    std::vector<MopsDetection> allDets;
    std::vector<Track> allTracks;
    unsigned int curIndex = 0;
    std::vector<std::vector< double > > positions;
    std::vector<std::vector<double> > answerKey;
    
    // build a vector of positions with 
    for (unsigned int i = 0; i < NUM_CALLS; i++ ) {
        Track newTrack; 
        double p0, v, acc;
        p0 = (i % 100) / 50.;
        v = (i % 13) / 20. - 7.;
        acc = (i % 7) / 50. - 3.;
        if (i % 3 == 0) {
            v *= -1;
        }
        if (i % 2 == 0) {
            acc *= -1;
        }
        std::vector<double> motion;
        
        motion.push_back(p0);
        motion.push_back(v);
        motion.push_back(acc);
        answerKey.push_back(motion);
        
        
        for (unsigned int j = 0; j < times.size(); j++) {
            addDetectionAt(times.at(j), projectLoc(times.at(j), p0, v, acc),
                           projectLoc(times.at(j), p0, v, acc) + 5,  allDets);
            newTrack.addDetection(curIndex, allDets);
            curIndex += 1;
        }
    }
    
    for (unsigned int i = 0; i < allTracks.size(); i++) {
        allTracks.at(i).calculateBestFitQuadratic(allDets);
        double epoch, ra0, raV, raAcc, dec0, decV, decAcc;
        allTracks.at(i).getBestFitQuadratic(epoch, ra0, raV, raAcc,
                                            dec0, decV, decAcc);
        BOOST_CHECK(Eq(epoch, 0));
        BOOST_CHECK(Eq(ra0,     answerKey.at(i).at(0)));
        BOOST_CHECK(Eq(raV,     answerKey.at(i).at(1)));
        BOOST_CHECK(Eq(raAcc,   answerKey.at(i).at(1)));
        BOOST_CHECK(Eq(dec0,     answerKey.at(i).at(0) + 5));
        BOOST_CHECK(Eq(decV,     answerKey.at(i).at(1)));
        BOOST_CHECK(Eq(decAcc,   answerKey.at(i).at(1)));

    }
}




}} // close lsst::mops

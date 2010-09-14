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

// for rand()
#include <cstdlib> 
// for printing timing info
#include <time.h>



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





    }} // close lsst::mops

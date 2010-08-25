#define BOOST_TEST_MODULE collapseTracklets

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
#include "lsst/mops/fileUtils.h"
#include "lsst/mops/rmsLineFit.h"
#include "lsst/mops/removeSubsets.h"
#include "lsst/mops/daymops/collapseTrackletsAndPostfilters/collapseTracklets.h"


bool Eq(double a, double b) 
{
    double epsilon = 1e-10;
    return (fabs(a - b) < epsilon);
}










using namespace lsst::mops;




///////////////////////////////////////////////////////////////////////
//    COLLAPSETRACKLETS TESTS
///////////////////////////////////////////////////////////////////////



BOOST_AUTO_TEST_CASE( isSane_blackbox_1 )
{
    std::vector<Tracklet> pairs;
    // first try calling with empty pairs
    BOOST_CHECK(isSane(100, &pairs) == true);

    //create some tracklets, put them in pairs
    for (int i = 0; i < 100; i++) {        
        Tracklet tmp;
        for (int j = 0; j < i; j++) {
            tmp.indices.insert(j);
        }
        pairs.push_back(tmp);
    }

    BOOST_CHECK(isSane(100, &pairs) == true);
    
    BOOST_CHECK(isSane(99, &pairs) == true);
        
}




// this function is now 'hidden' to slim down collapseTracklets.h.

// BOOST_AUTO_TEST_CASE( trackletsAreCompatible_blackbox_1 )
// {

//     std::vector<MopsDetection> dets;
//     MopsDetection tmpDet;
//     // 4 matching detections
//     tmpDet.fromMITIString("0 5330.0 1.0 1.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("0 5330.1 0.0 0.0 22.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("0 5330.2 359. 359. 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("0 5330.3 358. 358. 22.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);

//     // 4 unrelated detections
//     tmpDet.fromMITIString("0 5330.0 1.0 .8 20.0 1337 dummy2 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("0 5330.1 0.0 .2 22.0 1337 dummy2 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("0 5330.2 359. 358.8 20.0 1337 dummy2 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("0 5330.3 358. 358.2 22.0 1337 dummy2 0.0 0.0"); 
//     dets.push_back(tmpDet);

//     Tracklet t1, t2;
//     t1.indices.insert(0);
//     t1.indices.insert(1);  
//     t2.indices.insert(4);
//     t2.indices.insert(5);

//     BOOST_CHECK(trackletsAreCompatible(&dets, t1, t2) == false);
//     BOOST_CHECK(trackletsAreCompatible(&dets, t2, t1) == false);

//     t1.indices.clear();
//     t2.indices.clear();

//     t2.indices.insert(0);
//     t2.indices.insert(1);  
//     t1.indices.insert(4);
//     t1.indices.insert(5);
//     BOOST_CHECK(trackletsAreCompatible(&dets, t2, t1) == false);

//     t1.indices.clear();
//     t2.indices.clear();
 
//     t1.indices.insert(0);
//     t1.indices.insert(1);  
//     t1.indices.insert(2);  
//     t2.indices.insert(0);
//     t2.indices.insert(3);
//     BOOST_CHECK(trackletsAreCompatible(&dets, t1, t2) == true);
//     BOOST_CHECK(trackletsAreCompatible(&dets, t2, t1) == true);


//     t1.indices.clear();
//     t2.indices.clear();
    
//     t2.indices.insert(4);
//     t2.indices.insert(5);
//     BOOST_CHECK(trackletsAreCompatible(&dets, t1, t2) == true);
//     BOOST_CHECK(trackletsAreCompatible(&dets, t2, t1) == true);
    
// }






BOOST_AUTO_TEST_CASE( collapse_blackbox_2)
{
    Tracklet t1, t2;
    t1.indices.insert(1);
    t1.indices.insert(2);
    t1.indices.insert(3);
    t2.indices.insert(3);
    t2.indices.insert(4);

    TrackletCollapser myTC;
    myTC.collapse(t1, t2);
    BOOST_CHECK(t2.indices.size() == 4);
    BOOST_CHECK(t2.isCollapsed == true);
    BOOST_CHECK(t1.isCollapsed == true);
}





BOOST_AUTO_TEST_CASE( leastSquaresSolveForRADecLinear_blackbox_1)
{
    std::vector<MopsDetection> dets;
    MopsDetection tmpDet;

    tmpDet.fromMITIString("0 5330.0 10.0 10.0 20.0 1337 dummy 0.0 0.0"); 
    dets.push_back(tmpDet);
    tmpDet.fromMITIString("0 5331.0 11.0 11.0 20.0 1337 dummy 0.0 0.0");  
    dets.push_back(tmpDet);
    tmpDet.fromMITIString("0 5332.0 12.0 12.0 20.0 1337 dummy 0.0 0.0"); 
    dets.push_back(tmpDet);

    std::vector<double> RASlopeAndOffset, DecSlopeAndOffset;
    leastSquaresSolveForRADecLinear(&dets, RASlopeAndOffset, DecSlopeAndOffset,
						5330.0);
    BOOST_REQUIRE(RASlopeAndOffset.size() == 2);
    BOOST_REQUIRE(DecSlopeAndOffset.size() == 2);

    BOOST_CHECK(Eq(RASlopeAndOffset[0], 1.));
    BOOST_CHECK(Eq(RASlopeAndOffset[1], 10.));
    BOOST_CHECK(Eq(DecSlopeAndOffset[0], 1.));
    BOOST_CHECK(Eq(DecSlopeAndOffset[1], 10.));

}





// this function is also 'hidden' now...

// BOOST_AUTO_TEST_CASE( getAvgSqDist_blackbox_1)
// {
    
    
//     std::vector<MopsDetection> dets;
//     MopsDetection tmpDet;

//     tmpDet.fromMITIString("0 5330.0 10.0 10.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("0 5331.0 11.0 11.0 20.0 1337 dummy 0.0 0.0");  
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("0 5332.0 12.0 12.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);

//     std::vector<double> RASlopeAndOffset, DecSlopeAndOffset;
//     RASlopeAndOffset.push_back(1.0);
//     RASlopeAndOffset.push_back(10.0);
//     DecSlopeAndOffset.push_back(1.0);
//     DecSlopeAndOffset.push_back(10.0);

//     Tracklet t;
//     t.indices.insert(0);
//     t.indices.insert(1);
//     t.indices.insert(2);

//     double sqDist = getAverageSqDist(RASlopeAndOffset, DecSlopeAndOffset,
//                                                         &dets, &t, 5330.0);
//     BOOST_CHECK(Eq(sqDist, 0.0));
    

// }



// BOOST_AUTO_TEST_CASE( getAvgSqDist_blackbox_2)
// {
    
    
//     std::vector<MopsDetection> dets;
//     MopsDetection tmpDet;

//     // note that detections are each 1 degree off of the ideal
//     tmpDet.fromMITIString("0 5330.0 10.0 11.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("0 5331.0 11.0 10.0 20.0 1337 dummy 0.0 0.0");  
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("0 5332.0 12.0 13.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);

//     std::vector<double> RASlopeAndOffset, DecSlopeAndOffset;
//     RASlopeAndOffset.push_back(1.0);
//     RASlopeAndOffset.push_back(10.0);
//     DecSlopeAndOffset.push_back(1.0);
//     DecSlopeAndOffset.push_back(10.0);

//     Tracklet t;
//     t.indices.insert(0);
//     t.indices.insert(1);
//     t.indices.insert(2);

//     double sqDist = getAverageSqDist(RASlopeAndOffset, DecSlopeAndOffset,
//                                                         &dets, &t, 5330.0);
//     BOOST_CHECK(Eq(sqDist, 1.0));
    

// }


/* writeTrackletsToOutFile: leaving this as TBD for now.  interacting with files
   is generally (in my limited experience) not considered smart unit testing
   behavior...
*/








BOOST_AUTO_TEST_CASE( doCollapsingPopulateOutputVector_blackbox_1)
{

    std::vector<MopsDetection> dets;
    MopsDetection tmpDet;

    tmpDet.fromMITIString("0 5330.0 10.0 10.0 20.0 1337 dummy 0.0 0.0"); 
    dets.push_back(tmpDet);
    tmpDet.fromMITIString("1 5331.0 11.0 11.0 20.0 1337 dummy 0.0 0.0");  
    dets.push_back(tmpDet);
    tmpDet.fromMITIString("2 5332.0 12.0 12.0 20.0 1337 dummy 0.0 0.0"); 
    dets.push_back(tmpDet);
    tmpDet.fromMITIString("3 5333.0 13.0 13.0 20.0 1337 dummy 0.0 0.0"); 
    dets.push_back(tmpDet);

    std::vector<double> tolerances;
    tolerances.push_back(.01);
    tolerances.push_back(.01);
    tolerances.push_back(1.);
    tolerances.push_back(.01);

    std::vector<Tracklet> pairs;
    Tracklet tmpTracklet;
    tmpTracklet.indices.insert(0);
    tmpTracklet.indices.insert(1);
    pairs.push_back(tmpTracklet);
    tmpTracklet.indices.clear();
    tmpTracklet.indices.insert(2);
    tmpTracklet.indices.insert(3);
    pairs.push_back(tmpTracklet);

    std::vector<Tracklet> output;

    TrackletCollapser myTC;
    myTC.doCollapsingPopulateOutputVector(&dets, pairs, tolerances, output, false, false, false, 0.0, false);
    std::set <unsigned int> allIndices;
    allIndices.insert(0);
    allIndices.insert(1);
    allIndices.insert(2);
    allIndices.insert(3);
    BOOST_REQUIRE(output.size() == 1);
    BOOST_CHECK(output[0].indices == allIndices);
        
}




BOOST_AUTO_TEST_CASE( doCollapsingPopulateOutputVector_blackbox_2)
{

    std::vector<MopsDetection> dets;
    MopsDetection tmpDet;

    //slightly fuzzier data as well as fuzzier search
    tmpDet.fromMITIString("0 5330.0 10.0 10.0 20.0 1337 dummy 0.0 0.0"); 
    dets.push_back(tmpDet);
    tmpDet.fromMITIString("1 5331.0 11.0 11.5 20.0 1337 dummy 0.0 0.0");  
    dets.push_back(tmpDet);
    tmpDet.fromMITIString("2 5332.0 12.0 12.0 20.0 1337 dummy 0.0 0.0"); 
    dets.push_back(tmpDet);
    tmpDet.fromMITIString("3 5333.0 13.0 13.5 20.0 1337 dummy 0.0 0.0"); 
    dets.push_back(tmpDet);

    std::vector<double> tolerances;
    tolerances.push_back(2);
    tolerances.push_back(2);
    tolerances.push_back(45);
    tolerances.push_back(1);

    std::vector<Tracklet> pairs;
    Tracklet tmpTracklet;
    tmpTracklet.indices.insert(0);
    tmpTracklet.indices.insert(1);
    pairs.push_back(tmpTracklet);
    tmpTracklet.indices.clear();
    tmpTracklet.indices.insert(2);
    tmpTracklet.indices.insert(3);
    pairs.push_back(tmpTracklet);

    std::vector<Tracklet> output;

    TrackletCollapser myTC;
    myTC.doCollapsingPopulateOutputVector(&dets, pairs, tolerances, output, false, false, false, 0.0, false);
    std::set <unsigned int> allIndices;
    allIndices.insert(0);
    allIndices.insert(1);
    allIndices.insert(2);
    allIndices.insert(3);
    BOOST_REQUIRE(output.size() == 1);
    BOOST_CHECK(output[0].indices == allIndices);
        
}




BOOST_AUTO_TEST_CASE( doCollapsingPopulateOutputVector_blackbox_3)
{

    std::vector<MopsDetection> dets;
    MopsDetection tmpDet;

    //slightly fuzzier data, strict search - should not collapse anything.
    tmpDet.fromMITIString("0 5330.0 10.0 10.0 20.0 1337 dummy 0.0 0.0"); 
    dets.push_back(tmpDet);
    tmpDet.fromMITIString("1 5331.0 11.0 11.5 20.0 1337 dummy 0.0 0.0");  
    dets.push_back(tmpDet);
    tmpDet.fromMITIString("2 5332.0 12.0 12.0 20.0 1337 dummy 0.0 0.0"); 
    dets.push_back(tmpDet);
    tmpDet.fromMITIString("3 5333.0 13.0 13.5 20.0 1337 dummy 0.0 0.0"); 
    dets.push_back(tmpDet);

    std::vector<double> tolerances;
    tolerances.push_back(.5);
    tolerances.push_back(.5);
    tolerances.push_back(5);
    tolerances.push_back(.5);

    std::vector<Tracklet> pairs;
    Tracklet tmpTracklet;
    tmpTracklet.indices.insert(0);
    tmpTracklet.indices.insert(1);
    pairs.push_back(tmpTracklet);
    tmpTracklet.indices.clear();
    tmpTracklet.indices.insert(2);
    tmpTracklet.indices.insert(3);
    pairs.push_back(tmpTracklet);

    std::vector<Tracklet> output;

    TrackletCollapser myTC;
    myTC.doCollapsingPopulateOutputVector(&dets, pairs, tolerances, output, false, false, false, 0.0, false);
    std::set <unsigned int> allIndices;
    allIndices.insert(0);
    allIndices.insert(1);
    allIndices.insert(2);
    allIndices.insert(3);
    BOOST_CHECK(output.size() == 2);
        
}





// this function is now hidden.  It should probably be tested by whitebox testing of doCollapsing.
// BOOST_AUTO_TEST_CASE( detectionsHaveGT2UniqueTimes_blackbox_1)
// {
//     std::vector<MopsDetection> dets;
//     MopsDetection tmpDet;

//     tmpDet.fromMITIString("0 5330.0 10.0 10.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("1 5331.0 11.0 11.5 20.0 1337 dummy 0.0 0.0");  
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("2 5332.0 12.0 12.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("3 5333.0 13.0 13.5 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
    
//     BOOST_CHECK(detectionsHaveGT2UniqueTimes(&dets) == true);

//     dets.clear();
//     tmpDet.fromMITIString("0 5330.0 10.0 10.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("1 5330.0 11.0 11.5 20.0 1337 dummy 0.0 0.0");  
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("2 5331.0 12.0 12.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("3 5331.0 13.0 13.5 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);

//     BOOST_CHECK(detectionsHaveGT2UniqueTimes(&dets) == false);


// }




// BOOST_AUTO_TEST_CASE( setPhysicalParamsVector_blackbox_1)
// {
//     std::vector<MopsDetection> dets;
//     MopsDetection tmpDet;

//     tmpDet.fromMITIString("0 5330.0 10.0 10.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("1 5331.0 11.0 11.0 20.0 1337 dummy 0.0 0.0");  
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("2 5332.0 12.0 12.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("3 5333.0 13.0 13.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);     

//     std::vector<double> physicalParams(4);
//     setPhysicalParamsVector(&dets, physicalParams, 5331.5);
//     BOOST_CHECK(Eq(physicalParams[0], 11.5));
//     BOOST_CHECK(Eq(physicalParams[1], 11.5));
//     BOOST_CHECK(Eq(physicalParams[2], 45));
//     BOOST_CHECK(Eq(physicalParams[3], sqrt(2)));
// }



// BOOST_AUTO_TEST_CASE( setPhysicalParamsVector_blackbox_2)
// {
//     std::vector<MopsDetection> dets;
//     MopsDetection tmpDet;

//     tmpDet.fromMITIString("0 5330.0 358.0 358.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("1 5331.0 359.0 359.0 20.0 1337 dummy 0.0 0.0");  
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("2 5332.0 0.0 0.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);
//     tmpDet.fromMITIString("3 5333.0 1.0 1.0 20.0 1337 dummy 0.0 0.0"); 
//     dets.push_back(tmpDet);     

//     std::vector<double> physicalParams(4);
//     setPhysicalParamsVector(&dets, physicalParams, 5331.5);
//     BOOST_CHECK(Eq(physicalParams[0], 359.5));
//     BOOST_CHECK(Eq(physicalParams[1], 359.5));
//     BOOST_CHECK(Eq(physicalParams[2], 45));
//     BOOST_CHECK(Eq(physicalParams[3], sqrt(2)));
// }







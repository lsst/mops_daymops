// -*- LSST-C++ -*-
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <math.h>
#include <sstream>
#include <time.h>

#include "lsst/mops/KDTree.h"
#include "lsst/mops/common.h"
#include "lsst/mops/daymops/fieldProximity/fieldProximity.h"

#define uint unsigned int

#define LEAF_NODE_SIZE 8

namespace lsst { namespace mops {




void addPAV(const FieldProximityPoint &p, 
            uint objId,
            std::vector<PointAndValue<uint> > &allTrackPoints)
{
    PointAndValue<uint> tmpPAV;
    std::vector<double> tmpPt;
    tmpPt.push_back(p.getEpochMJD());
    tmpPt.push_back(convertToStandardDegrees(p.getRA()));
    tmpPt.push_back(convertToStandardDegrees(p.getDec()));
    tmpPAV.setPoint(tmpPt);
    tmpPAV.setValue(objId);

    allTrackPoints.push_back(tmpPAV);
}
        



void buildPavsForEphem(const std::vector<FieldProximityTrack> &tracks, 
                       std::vector<PointAndValue<uint> > &allTrackPoints)
{
                         
    // build KD-Tree PointsAndValues to put in the tree.  Each
    // ephemeris (time, ra, dec) is a point in the tree and the
    // value is the index of the track holding that ephemeris.

    for(unsigned int i=0; i < tracks.size(); i++){
        const std::vector<FieldProximityPoint> * pointsHere;
        pointsHere = tracks.at(i).getPoints();
        for (unsigned int j=0; j < pointsHere->size(); j++)
        {
            addPAV(pointsHere->at(j), i, allTrackPoints);
        }
    }

}





/* Helper function for getProximity.  When a tree search reveals that
 * a given track may pass through a given image, we need to see the
 * points nearest that image and check whether they do fall within the
 * image itself.

 * modify args to point at the two points in the track which happen
 * immediately before and after the specified image, or return an
 * exception.  Note that this should be fine for our current data sets
 * (may 23 2011) since we pre-filter and make sure that any object
 * remotely near the images is included in the data set.
 *
 * this could be done more efficiently if we assumed that the track's
 * associated points were in sorted order.
 */
void getNearestPoints(const Field &img, 
                      const FieldProximityTrack &track,
                      FieldProximityPoint &before,
                      FieldProximityPoint &after)
{
    double imgMjd = img.getEpochMJD();

    const std::vector<FieldProximityPoint> * trackPts = track.getPoints();
    if (trackPts->size() < 2) {
        throw LSST_EXCEPT(BadParameterException, 
       "Field Proximity track contains < 2 points. This is impossible to work with.");
    }
    uint beforeIndex = 0;
    double maxBefore, minAfter;
    bool foundBefore = false;
    bool foundAfter = false;
    uint afterIndex = 0;
    for (uint i = 0; i < trackPts->size(); i++) {
        double mjdHere = trackPts->at(i).getEpochMJD();
        if (mjdHere < imgMjd) {
            if ((foundBefore == false) || 
                (maxBefore < mjdHere)) {
                foundBefore = true;
                maxBefore = mjdHere;
                beforeIndex = i;
            }
        }
        else if (mjdHere > imgMjd) {
            if ((foundAfter == false) || 
                (minAfter > mjdHere)) {
                foundAfter = true;
                minAfter = mjdHere;
                afterIndex = i;
            }
        }
    }
    if ((!foundBefore)||(!foundAfter)) {
        std::cerr << "FAILURE! Image mjd, ra, dec were " <<
            imgMjd << ", " << img.getRA() << ", " << img.getDec() << "\n";
        std::cerr << " Track for object " << track.getID() << " had times/locs: \n";
        for (uint i = 0; i < trackPts->size(); i++) {
            std::cerr << "\t" << trackPts->at(i).getEpochMJD()
                      << ", " << trackPts->at(i).getRA() 
                      << ", " << trackPts->at(i).getDec() << "\n";
        }
        throw LSST_EXCEPT(BadParameterException,
                          "Track did not have points before/after the image time.");
    }

    before = (trackPts->at(beforeIndex));
    after = (trackPts->at(afterIndex));
}



bool isInsideImage(const Field &img, const FieldProximityTrack &t) 
{
    FieldProximityPoint before, after;
    getNearestPoints(img, t, before, after);
    double ra0 = before.getRA();
    double dec0 = before.getDec();
    double ra1 = after.getRA();
    double dec1 = after.getDec();
    while (ra0 - ra1 > 180.) {
        ra1 += 360;
    }
    while (ra0 - ra1 < -180.) {
        ra1 -= 360.;
    }
    while (dec0 - dec1 > 180.) {
        dec1 += 360;
    }
    while (dec0 - dec1 < -180.) {
        dec1 -= 360.;
    }

    double t0 = before.getEpochMJD();
    double dt = after.getEpochMJD() - t0;
    
    double raSlope = (ra1 - ra0) / dt;
    double decSlope = (dec1 - dec0) / dt;
    
    double imgT = img.getEpochMJD();
    double predRa = ra0 + raSlope*(imgT - t0);
    double predDec = dec0 + decSlope*(imgT - t0);
    
    double dist = angularDistanceRADec_deg(predRa, predDec, img.getRA(), img.getDec());

    return (dist <= img.getRadius());
}



void getProximity(std::vector<std::pair<uint, uint> > &resultsVec,  
                  KDTree<uint> &myTree,
                  const std::vector<Field> &queryPoints,
                  const std::vector<FieldProximityTrack> &allTracks)
{
    // vectors of hyperRectangleSearch parameters
    std::vector<double> tolerances;
  
    std::vector<GeometryType> ephemGeometry;
    ephemGeometry.push_back(EUCLIDEAN);
    ephemGeometry.push_back(CIRCULAR_DEGREES);
    ephemGeometry.push_back(CIRCULAR_DEGREES);

    // hyperRectangleSearch result container
    std::vector<PointAndValue<uint> > queryResults;
    
    for (uint i = 0; i < queryPoints.size(); i++) {
        /* wonderfully confusing KDTree interface.  basically, find
         * any ephemeris which is inside the image and within 1 day of
         * the image.
         */

        std::vector<double> queryPoint(3,0);
        queryPoint[0] = queryPoints[i].getEpochMJD();
        queryPoint[1] = queryPoints[i].getRA();
        queryPoint[2] = queryPoints[i].getDec();
        std::vector<double> queryRange(3,0);
        queryRange[0] = 1;
        queryRange[1] = queryPoints[i].getRadius();
        queryRange[2] = queryPoints[i].getRadius();

        queryResults = myTree.hyperRectangleSearch(queryPoint, 
                                                   queryRange,
                                                   ephemGeometry);

        // we now know all the tracks which pass near this image.

        // we may get some redundant entries from the tree search;
        // avoid needless processing by removing redundant entries
        std::set<uint> queryResultsSet;
        for (uint j = 0; j < queryResults.size(); j++) {
            queryResultsSet.insert(queryResults[j].getValue());
        }
        std::set<uint>::const_iterator resultIter;
        for (resultIter = queryResultsSet.begin(); 
             resultIter != queryResultsSet.end(); 
             resultIter++) {

            const FieldProximityTrack* matchingTrack = 
                &(allTracks.at(*resultIter));
            if (isInsideImage(queryPoints[i], *matchingTrack)) {
                resultsVec.push_back(std::make_pair(queryPoints[i].getFieldID(),
                                                    matchingTrack->getID()));

            }
        }
    }

}





bool compByTime(const Field &f1, const Field &f2)
{
    return f1.getEpochMJD() < f2.getEpochMJD();
}

void fieldProximity(const std::vector<FieldProximityTrack> &allTracks,
                    std::vector<Field> &queryFields,
                    std::vector<std::pair<unsigned int, unsigned int>  > &results,
                    double distThresh)
{

    if(queryFields.size() > 0 && allTracks.size() > 0){
        
        std::sort(queryFields.begin(), queryFields.end(), compByTime);
        
        std::vector<PointAndValue<uint> > allEphem;
        time_t currentTime;
        time(&currentTime);
        std::cout << "Massaging data for tree construction at " 
                  << ctime(&currentTime) << "\n";
        buildPavsForEphem(allTracks, allEphem);
        // build a tree of (time, ra, dec) -> track index 
        time(&currentTime);
        std::cout << "Building tree at " 
                  << ctime(&currentTime) << "\n";
        KDTree<uint> myTree(allEphem, 3, LEAF_NODE_SIZE);        
        time(&currentTime);
        std::cout << "Searching tree at " 
                  << ctime(&currentTime) << "\n";
        getProximity(results, myTree, queryFields, allTracks);
        time(&currentTime);
        std::cout << "Finished searching at " 
                  << ctime(&currentTime) << "\n";
        
    }
    
}


// legacy interface
std::vector<std::pair<unsigned int, unsigned int>  > 
fieldProximity(const std::vector<FieldProximityTrack> &allTracks,
               std::vector<Field> &queryFields,
               double distThresh)
{
    std::vector<std::pair<unsigned int, unsigned int>  > toRet;
    fieldProximity(allTracks, queryFields, toRet, distThresh);
    return toRet;
}











}} // close lsst::mops 

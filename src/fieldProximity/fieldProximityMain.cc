// -*- LSST-C++ -*-
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <getopt.h>
#include <map>
#include <string>
#include <sstream>
#include <utility>
#include <sstream>

#include "lsst/mops/daymops/fieldProximity/fieldProximity.h"
#include "lsst/mops/daymops/fieldProximity/TrackForFieldProximity.h"

#define uint unsigned int

    
namespace lsst { namespace mops {

void debugPrintTracks(const std::vector<FieldProximityTrack> &allTracks)
{
    std::vector<FieldProximityTrack>::const_iterator i;
    for (i = allTracks.begin(); i != allTracks.end(); i++) {
        std::cout << "Track " << i->getID() << " : \n";

        const std::vector<FieldProximityPoint> * points = i->getPoints();
        std::vector<FieldProximityPoint>::const_iterator j;
        
        for (j = points->begin(); j != points->end(); j++) {
            std::cout << std::setprecision(10);
            std::cout << "\tmjd=" << j->getEpochMJD() << 
                " pos=(" << j->getRA() << ", " << j->getDec() << ")\n";
        }
    }
}

void debugPrintFields(const std::vector<Field> &allFields)
{
    std::vector<Field>::const_iterator i;
    for (i = allFields.begin(); i != allFields.end(); i++) {
        std::cout << std::setprecision(10);
        std::cout << "Field " << i->getFieldID() << " : mjd=" << 
            i->getEpochMJD() << " pos=(" << i->getRA() <<
            ", " << i->getDec() << ") radius=" << i->getRadius() << "\n";
                  
    }
}


void populateFieldsVec(std::string fieldsFileName, std::vector<Field> &queryFields)
{
    std::string line;
    std::ifstream fieldsFile;
    fieldsFile.open(fieldsFileName.c_str());
    if (!fieldsFile.is_open()) {
         throw LSST_EXCEPT(FileException,
                           "Failed to open fields file " 
                           + fieldsFileName + " - does this file exist?\n");
     }

    std::getline(fieldsFile, line);
    while (fieldsFile.fail() == false) {
         uint tmpFieldId;
         double tmpExpMjd;
         double tmpRa;
         double tmpDec;
         double tmpRadius;
         std::istringstream ss(line);
         ss.exceptions(std::ifstream::failbit | std::ifstream::badbit);    
         ss >> tmpFieldId;
         ss >> std::ws;
         ss >> tmpExpMjd;
         ss >> std::ws;
         ss >> tmpRa;
         ss >> std::ws;
         ss >> tmpDec;
         ss >> std::ws;
         ss >> tmpRadius;
         ss >> std::ws;

         Field tmpField;
         tmpField.setFieldID(tmpFieldId);
         tmpField.setEpochMJD(tmpExpMjd);
         tmpField.setRA(tmpRa);
         tmpField.setDec(tmpDec);
         tmpField.setRadius(tmpRadius);

         queryFields.push_back(tmpField);
         std::getline(fieldsFile, line);
     }
    
}



/* helper function for buildTrackVector. Add the point to the
 * right track; adding entries in the map/vector as needed. */
void addPoint(uint trackId, 
              FieldProximityPoint toAdd, 
              std::vector<FieldProximityTrack> &allTracks,
              std::map<uint, uint> &trackMap)
{
    if (trackMap.find(trackId) == trackMap.end()) {
        // create a new track
        FieldProximityTrack newTrack;
        newTrack.setID(trackId);
        newTrack.addPoint(toAdd);
        allTracks.push_back(newTrack);
        trackMap[trackId] = allTracks.size() - 1;
    }
    else {
        // just update this entry
        allTracks.at(trackMap[trackId]).addPoint(toAdd);
    }
}



void buildTrackVector(std::string tracksFileName, 
                      std::vector<FieldProximityTrack> &allTracks)
{
    /* trackMap will hold a mapping from track ID to index in the
     * allTracks vector.  This is needed because not all points in a
     * track are contiguous in the file, so we need to do lookups. */
    std::map<uint, uint> trackMap; 
    std::string line;
    std::ifstream tracksFile;
    tracksFile.open(tracksFileName.c_str());
     if (!tracksFile.is_open()) {
         throw LSST_EXCEPT(FileException,
                           "Failed to open tracks file " 
                           + tracksFileName + " - does this file exist?\n");
     }

     std::getline(tracksFile, line);
     while (tracksFile.fail() == false) {
         std::istringstream ss(line);
         ss.exceptions(std::ifstream::failbit | std::ifstream::badbit);    
         unsigned int tmpId;
         double tmpMjd;
         double tmpRa;
         double tmpDec;
         FieldProximityPoint tmpPoint;
         ss >> tmpId;
         ss >> std::ws;
         ss >> tmpMjd;
         ss >> std::ws;
         ss >> tmpRa;
         ss >> std::ws;
         ss >> tmpDec;
         ss >> std::ws;
         tmpPoint.setRA(tmpRa);
         tmpPoint.setDec(tmpDec);
         tmpPoint.setEpochMJD(tmpMjd);
         
         addPoint(tmpId, tmpPoint, allTracks, trackMap);
         std::getline(tracksFile, line);
     }
}

void writeFieldMatches(
    std::string outFileName,
    const std::vector<std::pair<unsigned int, unsigned int> > &matches, 
    const std::vector<Field> &allFields, 
    const std::vector<FieldProximityTrack> &allTracks)
{
    std::ofstream outFile;
    outFile.open(outFileName.c_str());
    if (!outFile.is_open()) {
        throw LSST_EXCEPT(FileException,
     "Failed to open output file " + outFileName + " - do you have permission?\n");
     }
    
    for (uint i = 0; i < matches.size(); i++) {
        std::pair<uint, uint> match = matches.at(i);
        outFile << match.first << 
            "\t" << match.second << "\n";;
        
    }


}


}} // close lsst::mops




int main(int argc, char *args[])
{

    std::string fieldsFile, tracksFile, outFile = "";
    double maxDist = .01;
    
    std::string USAGE = "Usage: fieldProximity -f <fields file> -t <tracks file> -o <output file> [-r <threshold degrees>]\n";
  
    if(argc < 3){
        std::cout << USAGE;
        exit(1);
    }

    static const struct option longOpts[] = {
        { "fieldsFile", required_argument, NULL, 'f' },
        { "tracksFile", required_argument, NULL, 't' },
        { "outFile", required_argument, NULL, 'o' },
        { "distThresh", required_argument, NULL, 'r' },
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
    };


    int longIndex = -1;
    const char *optString = "f:t:o:r:h";
    int opt = getopt_long( argc, args, optString, longOpts, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
        case 'f':
            fieldsFile = std::string(optarg); 
            break;
        case 't':
            tracksFile = std::string(optarg);
            break;
        case 'o':
            outFile = std::string(optarg);
            break;
        case 'r':
            maxDist = atof(optarg);
            break;
        case 'h':
            std::cout << USAGE;
            exit(0);
        default:
            break;
        }
        opt = getopt_long( argc, args, optString, longOpts, &longIndex );
    }

    if ((fieldsFile == "") || (tracksFile == "") || (outFile == "")) {
        std::cerr << USAGE;
        exit(1);
    }

    //vector of RA and dec pairs for later searching
    std::vector<lsst::mops::Field> queryFields; 

    std::vector<lsst::mops::FieldProximityTrack > allTracks;

    std::cout << "Reading list of fields from " << fieldsFile << "\n";
    lsst::mops::populateFieldsVec(fieldsFile, queryFields);

    std::cout << "Reading list of ephemerides from " << tracksFile << "\n";
    lsst::mops::buildTrackVector(tracksFile, allTracks);

    //debugPrintFields(queryFields);
    //debugPrintTracks(allTracks);

    std::vector<std::pair<unsigned int, unsigned int> > matches;

    std::cout << "Looking for possible overlaps..."  << "\n";
    lsst::mops::fieldProximity(allTracks, queryFields, matches, maxDist);

    std::cout << "Writing results to " << outFile << "\n";
    lsst::mops::writeFieldMatches(outFile, matches, queryFields, allTracks);

    return 0;
}


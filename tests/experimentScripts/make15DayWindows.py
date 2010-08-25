#!/usr/bin/env python

""" jmyers

May 13 2010

A quickie script to build .dets files and .ids files for
linkTracklets, using sets of nightly dets and ids files.

Input dets should be in MITI format.

Input ids files should list detections as *DIASOURCE IDS*

Input and output files should have names like

1234.miti and 1234.maxv5.6.tracklets for the Dets file and Ids file for night 1234.



Output dets file will be in MITI format.

Output ids file will list detections *INDEXES INTO THE DETS FILE*
suitable for running C++ linkTracklets

"""

INPUT_DETS_DIR="/workspace1/jmyers/nightlyDiasAstromErr/"
INPUT_DETS_SUFFIX=".miti"
INPUT_TRACKLETS_DIR="/workspace1/jmyers/nightlyDiasAstromErr/tracklets/"
INPUT_TRACKLETS_SUFFIX=".maxv0.5.tracklets"
OUTPUT_DIR="/workspace1/jmyers/nightlyDiasAstromErr/15DayWindowsMaxv0.5/"
WINDOW_SIZE=15

import glob, sys


def addDetToIndexLookupTable(table, strDet, index):
    detId = int(strDet.split()[0])
    table[detId] = index



if __name__=="__main__":

    mitiFiles = glob.glob(INPUT_DETS_DIR + '*' + INPUT_DETS_SUFFIX)

    nights = sorted(map(lambda x: int(x[len(INPUT_DETS_DIR):-(len(INPUT_DETS_SUFFIX))]), mitiFiles))
    print "All nights: ", nights

    previousNights = []
    for night in nights:

        compatibleNights = sorted(filter(lambda x: x >= night and x <= night + WINDOW_SIZE, nights))

        #don't create rundundant work by creating identical sets
        if compatibleNights != previousNights:

            outDets = file(OUTPUT_DIR + "night_" + str(night) + "_through_" + str(compatibleNights[-1]) + ".dets", 'w')
            outIds =  file(OUTPUT_DIR + "night_" + str(night) + "_through_" + str(compatibleNights[-1]) + ".ids", 'w')

            outputDetsSize = 0
            

            for compatibleNight in compatibleNights:
                print "compatibleNights for night ", night, " :  ", compatibleNights
                print "adding ", compatibleNight, " to output set for ", night

                # add this night's detections and revised tracklets into outfile

                #use this to track det Id -> index into our new output file
                idToIndexTable = {}

                compatibleDets = file(INPUT_DETS_DIR + str(compatibleNight) + INPUT_DETS_SUFFIX, 'r')
                compatibleIds =  file(INPUT_TRACKLETS_DIR  + str(compatibleNight) + INPUT_TRACKLETS_SUFFIX, 'r')
                det = compatibleDets.readline()
                while det != "":
                    outDets.write(det)
                    addDetToIndexLookupTable(idToIndexTable, det, outputDetsSize)
                    outputDetsSize += 1
                    det = compatibleDets.readline()
                
                pair = compatibleIds.readline()
                while pair != "":
                    detIds = map(int, pair.split())
                    indexes = map(lambda x: idToIndexTable[x], detIds)
                    indexesStr = " ".join(map(str, indexes))
                    #print "writing output ", indexesStr
                    outIds.write(indexesStr + "\n")
                    pair = compatibleIds.readline()
                

            previousNights = compatibleNights
        
            

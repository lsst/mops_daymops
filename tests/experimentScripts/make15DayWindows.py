#!/bin/env python

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

INPUT_DETS_SUFFIX=".miti"
INPUT_TRACKLETS_SUFFIX=".maxv0.5.tracklets"
WINDOW_SIZE=15
# can specify these next 3 values on the command line instead
INPUT_DETS_DIR="/workspace1/jmyers/nightlyDiasAstromErr/"
INPUT_TRACKLETS_DIR="/workspace1/jmyers/nightlyDiasAstromErr/tracklets/"
OUTPUT_DIR="/workspace1/jmyers/nightlyDiasAstromErr/15DayWindowsMaxv0.5/"

# use new style 'windows' where 'nightNum' of window is end of data
# and extends BACKWARDS in time WINDOW_SIZE nights (instead of forward). 
NEW_WINDOWS = True

import glob, sys, os


def addDetToIndexLookupTable(table, strDet, index):
    detId = int(strDet.split()[0])
    table[detId] = index



if __name__=="__main__":

    if len(sys.argv)>2:
        INPUT_DETS_DIR = sys.argv[1]
        if INPUT_DETS_DIR != ".":
            INPUT_DETS_DIR = INPUT_DETS_DIR + "/"
        INPUT_TRACKLETS_DIR = sys.argv[2]
        if INPUT_TRACKLETS_DIR !=  ".":
            INPUT_TRACKLETS_DIR = INPUT_TRACKLETS_DIR + "/"
        OUTPUT_DIR = sys.argv[3]
        if OUTPUT_DIR != ".":
            OUTPUT_DIR = OUTPUT_DIR + "/"
        print "Reading MITI files from %s" %(INPUT_DETS_DIR)
        print "Reading tracklet files from %s" %(INPUT_TRACKLETS_DIR)
        print "Will write 15 day window *.dets and *.ids files to %s" %(OUTPUT_DIR)

    # Find the list of all *miti (OR INPUT_DETS_SUFFIX) files in the INPUT_DETS_DIR.
    mitiFiles = glob.glob(os.path.join(INPUT_DETS_DIR + '*' + INPUT_DETS_SUFFIX))

    # Sort the files in MJD order. 
    nights = sorted(map(lambda x: int(x[len(INPUT_DETS_DIR):-(len(INPUT_DETS_SUFFIX))]), mitiFiles))


    print "All nights: ", nights

    
    previousNights = []
    for night in nights:
        # Find nights which are within WINDOW_SIZE nights previous to 'night'. 
        if NEW_WINDOWS:
            compatibleNights = sorted(filter(lambda x: x <= night 
                                             and x >= (night - WINDOW_SIZE), nights))
        else:
            compatibleNights = sorted(filter(lambda x: x >= night 
                                             and x <= (night + WINDOW_SIZE), nights))
        print compatibleNights
        # Don't create redundant work by creating identical sets, and don't create output if
        # linkTracklets won't run (linkTracklets requires 3 nights of tracklets). 
        if (compatibleNights != previousNights) & (len(compatibleNights)>=3) : 
            if NEW_WINDOWS:                
                outDets = file(OUTPUT_DIR +
                               "linkTracklets_input_" + str(compatibleNights[-1]) + ".dets", "w")
                outIds = file(OUTPUT_DIR +
                              "linkTracklets_input_" + str(compatibleNights[-1]) + ".ids", "w")
            else:
                outDets = file(OUTPUT_DIR + "night_" + str(night) + 
                               "_through_" + str(compatibleNights[-1]) + ".dets", 'w')
                outIds =  file(OUTPUT_DIR + "night_" + str(night) + 
                               "_through_" + str(compatibleNights[-1]) + ".ids", 'w')

            outputDetsSize = 0
            
            for compatibleNight in compatibleNights:
                print "compatibleNights for night ", night, " :  ", compatibleNights
                print "adding ", compatibleNight, " to linkTracklets input set for ", night

                # Add this night's detections and revised tracklets into outfile.
                # Use this to track det Id -> index into our new output file. 
                idToIndexTable = {}

                # Open the tracklet input and output files.
                compatibleDets = file(INPUT_DETS_DIR + str(compatibleNight) + 
                                      INPUT_DETS_SUFFIX, 'r')
                compatibleIds =  file(INPUT_TRACKLETS_DIR + str(compatibleNight) + 
                                      INPUT_TRACKLETS_SUFFIX, 'r')
                # Read each diasource input into findTracklets. 
                det = compatibleDets.readline()
                while det != "":
                    outDets.write(det)
                    addDetToIndexLookupTable(idToIndexTable, det, outputDetsSize)
                    outputDetsSize += 1
                    det = compatibleDets.readline()
                # Read each tracklet found from findTracklets. 
                pair = compatibleIds.readline()
                while pair != "":
                    detIds = map(int, pair.split())
                    indexes = map(lambda x: idToIndexTable[x], detIds)
                    indexesStr = " ".join(map(str, indexes))
                    #print "writing output ", indexesStr
                    outIds.write(indexesStr + "\n")
                    pair = compatibleIds.readline()
                

            previousNights = compatibleNights
        
            
# Yes, this could be done more efficiently. 
# However, the final 'operations' will be fairly different and this is easy to code.
# Remember, this is writing the indexes into each file (line numbers of diasources in miti output)
# into the input files for linkTracklets. 

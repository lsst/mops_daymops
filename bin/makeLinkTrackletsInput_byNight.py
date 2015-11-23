#!/usr/bin/env python

""" jmyers

May 13 2010

A quickie script to build .dets files and .ids files for
linkTracklets, using sets of nightly dets and ids files.

Input ids files should list detections as *DIASOURCE IDS*

Input and output files should have names like

1234.dias and 1234.tracklets.final.byDiaIds for the Dets file and Ids file for night 1234.

TBGD: use the .final.byIndices format and you could cleverly speed this script up a lot.


Output ids file will list detections *INDEXES INTO THE DETS FILE*
suitable for running C++ linkTracklets

"""

INPUT_DETS_SUFFIX=".dias"
INPUT_TRACKLETS_SUFFIX=".tracklets.final.byDiaIds"
# can specify these next 3 values on the command line instead
INPUT_DETS_DIR="../../"
INPUT_TRACKLETS_DIR="../"
OUTPUT_DIR="./"

# use new style 'windows' where 'nightNum' of window is end of data
# and extends BACKWARDS in time WINDOW_SIZE nights (instead of forward).
NEW_WINDOWS = False #jmyers: i messed with this script 4-19-2012 and probably broke this.

import glob, sys, os


def addDetToIndexLookupTable(table, strDet, index):
    detId = int(strDet.split()[0])
    table[detId] = index


def writeDetsAndIds(detsFile, idsFile, compatibleNights):
    outputDetsSize = 0

    for compatibleNight in compatibleNights:

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

    outDets.close()
    outIds.close()
    print "done with ", compatibleNights
    print ""

if __name__=="__main__":

    WINDOW_SIZE=int(sys.argv[1])

    if len(sys.argv) > 2:
        INPUT_DETS_DIR=sys.argv[2]
        INPUT_TRACKLETS_DIR=sys.argv[3]
        OUTPUT_DIR=sys.argv[4]

    print "Reading MITI files from %s" %(INPUT_DETS_DIR)
    print "Reading tracklet files from %s" %(INPUT_TRACKLETS_DIR)
    print "Will write %d day window *.dets and *.ids files to %s" %(WINDOW_SIZE, OUTPUT_DIR)

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
        # Don't create redundant work by creating identical sets, and don't create output if
        # linkTracklets won't run (linkTracklets requires 3 nights of tracklets).
        if not ((compatibleNights != previousNights) & (len(compatibleNights)>=3)) :
            print "No track could start on night ", night
            print "Compatible nights are: ", compatibleNights
        else:
            # we can indeed do useful searching here.
            if NEW_WINDOWS:
                outDets = file(OUTPUT_DIR +
                               "linkTracklets_input_" + str(compatibleNights[-1]) + ".dets", "w")
                outIds = file(OUTPUT_DIR +
                              "linkTracklets_input_" + str(compatibleNights[-1]) + ".ids", "w")
            else:
                prefix = OUTPUT_DIR + "night_" + str(night) + "_through_" + str(compatibleNights[-1])
                outDets = file(prefix + ".dets", 'w')
                outIds =  file(prefix + ".ids", 'w')
                endTimeRange = compatibleNights[1]
                print "compatibleNights for night ", night, " :  ", compatibleNights
                print "adding ", compatibleNights, " to linkTracklets input set for ", night
                print "stop looking after finding all tracks starting before ", endTimeRange
                # jmyers: new: add start_t_range argument too
                tRangeFile = file(prefix + ".date.start_t_range", 'w')
                # jmyers: if using NEW_WINDOWS then probably set this to max?
                # and you need an end_t_range. so you'll have to fix
                # scripts elsewhere, too...
                tRangeFile.write(str(float(endTimeRange)))
                tRangeFile.close()

            writeDetsAndIds(outDets, outIds, compatibleNights)

            previousNights = compatibleNights



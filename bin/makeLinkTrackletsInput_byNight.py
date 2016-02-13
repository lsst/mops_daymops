#!/usr/bin/env python

""" jmyers

May 13 2010

A quickie script to build .dets files and .ids files for
linkTracklets, using sets of nightly dets and ids files.

Input ids files should list detections as *DIASOURCE IDS*

Input and output files should have names like

1234.dias and 1234.tracklets.final.byDiaIds for the Dets file and Ids file for night 1234.

TBGD: use the .final.byIndices format and you could cleverly speed this script up a lot.

February 2016 [moeyensj]: 
added in argparse including parameter specification, removed broken backwards extended window format

Output ids file will list detections *INDEXES INTO THE DETS FILE*
suitable for running C++ linkTracklets

"""
# can be specified on command line
INPUT_DETS_SUFFIX=".dias"
INPUT_TRACKLETS_SUFFIX=".tracklets.final.byDiaIds"
INPUT_DETS_DIR="../../"
INPUT_TRACKLETS_DIR="../"
OUTPUT_DIR="./"
WINDOW_SIZE = 15

import glob
import sys
import os

def addDetToIndexLookupTable(table, strDet, index):
    detId = int(strDet.split()[0])
    table[detId] = index

def writeDetsAndIds(outDets, outIds, compatibleNights, diasDir=INPUT_DETS_DIR, trackletDir=INPUT_TRACKLETS_DIR, diasSuffix=INPUT_DETS_SUFFIX, trackletSuffix=INPUT_TRACKLETS_SUFFIX):
    outputDetsSize = 0

    for compatibleNight in compatibleNights:

        # Add this night's detections and revised tracklets into outfile.
        # Use this to track det Id -> index into our new output file.
        idToIndexTable = {}

        # Open the tracklet input and output files.
        compatibleDets = open(diasDir + str(compatibleNight) +
                              diasSuffix, 'r')
        compatibleIds = open(trackletDir + str(compatibleNight) +
                              trackletSuffix, 'r')
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
            outIds.write(indexesStr + " \n")
            pair = compatibleIds.readline()

    outDets.close()
    outIds.close()
    print "done with ", compatibleNights
    print ""

if __name__=="__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Combines nightly DiaSource files into detection and tracklet files covering a window of nights.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("diasDir", type=str, help="Directory containing DiaSources (detections) seperated in files by night.")
    parser.add_argument("trackletDir", type=str,
                        help="Directory containing the tracklet files to be used for linking.")
    parser.add_argument("outputDir", type=str,
                        help="Directory to put output detection and ids files in.")
    parser.add_argument("--windowSize", type=int, default=WINDOW_SIZE,
                        help="The number of nights in a window.")
    parser.add_argument("--diasSuffix", type=str, default=INPUT_DETS_SUFFIX,
                        help="Dias file suffix.")
    parser.add_argument("--trackletSuffix", type=str, default=INPUT_TRACKLETS_SUFFIX,
                        help="Tracklet file suffix.")

    args = parser.parse_args()

    print "Reading dias files from %s" % (args.diasDir)
    print "Reading tracklet files from %s" % (args.trackletDir)
    print "Will write %d day window *.dets and *.ids files to %s" % (args.windowSize , args.outputDir)

    # Find the list of all detection files
    diasFiles = glob.glob(os.path.join(args.diasDir + '*' + args.diasSuffix))

    # Sort the files in MJD order.
    nights = sorted(map(lambda x: int(x[len(args.diasDir):-(len(args.diasSuffix))]), diasFiles))

    print "All nights: ", nights

    previousNights = []
    for night in nights:
        # Don't create redundant work by creating identical sets, and don't create output if
        # linkTracklets won't run (linkTracklets requires 3 nights of tracklets).
        compatibleNights = sorted(x for x in nights if x >= night and x <= (night + args.windowSize))

        if not ((compatibleNights != previousNights) & (len(compatibleNights)>=3)) :
            print "No track could start on night ", night
            print "Compatible nights are: ", compatibleNights
        else:
            # we can indeed do useful searching here.
            prefix = args.outputDir + "night_" + str(night) + "_through_" + str(compatibleNights[-1])
            outDets = open(prefix + ".dets", 'w')
            outIds =  open(prefix + ".ids", 'w')
            endTimeRange = compatibleNights[1]
            print "compatibleNights for night ", night, " :  ", compatibleNights
            print "adding ", compatibleNights, " to linkTracklets input set for ", night
            print "stop looking after finding all tracks starting before ", endTimeRange
            # jmyers: new: add start_t_range argument too
            tRangeFile = file(prefix + ".date.start_t_range", 'w')
            tRangeFile.write(str(float(endTimeRange)))
            tRangeFile.close()

        writeDetsAndIds(outDets, outIds, compatibleNights, diasDir=args.diasDir, trackletDir=args.trackletDir, diasSuffix=args.diasSuffix, trackletSuffix=args.trackletSuffix)

        previousNights = compatibleNights

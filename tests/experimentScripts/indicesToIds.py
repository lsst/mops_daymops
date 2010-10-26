#!/usr/bin/env python

""" jmyers july 15 2010

This script is for those who want to know what DiaSources go with which Track.

the current version of C linkTracklets "loses" the IDs of the
DiaSources used as input when we construct .tracklet.miti files from
the tracklets and detections.

Then, our new version of C LinkTracklets expresses its tracks as sets
of line numbers (array indices) into the .tracklet.miti files.  

use ../../makeCheatSheetForMITIFiles.py to get the Ids and MJDs of
DiaSources used at each of the lines in the .tracklet.miti files.
Then use this script to convert from the output of C linkTracklets
(sets of line numbers) over to sets of DiaSourceIds.

"""


DIA_IDS_DIR="/workspace0/jmyers/nightlyDiasAstromErr_linkTrackletsInfiles_maxv0.5_mint_0.01_30dayWindows/"
DIA_IDS_SUFFIX=".miti.diaIds"


TRACKS_DIR="/workspace0/jmyers/nightlyDiasAstromErr_linkTrackletsInfiles_maxv0.5_mint_0.01_30dayWindows/results/finished/"
TRACKS_SUFFIX=".c.tracks.byIndices"

import glob

TRACKFILES=glob.glob(TRACKS_DIR + "*" + TRACKS_SUFFIX)

TRACK_BY_DIA_SUFFIX=".c.tracks.byDiaId"



def lineNumToIdsAndMjds(detsFile):
    line = detsFile.readline()
    ids = []
    #mjds = []
    while line != "":
        diaId = int(line.split()[0])
        #mjd = float(line.split()[1])
        ids.append(diaId)        
        #mjds.append(mjd)
        line = detsFile.readline()
    return ids #, mjds











if __name__ == "__main__":

    import sys


    trackId = 0


    for trackFileName in TRACKFILES:

        diaIdsFileName = DIA_IDS_DIR + trackFileName[len(TRACKS_DIR):-len(TRACKS_SUFFIX)] + DIA_IDS_SUFFIX
        outFileName = TRACKS_DIR + trackFileName[len(TRACKS_DIR):-len(TRACKS_SUFFIX)] + TRACK_BY_DIA_SUFFIX

        print trackFileName, diaIdsFileName, outFileName


        trackFile = file(trackFileName, 'r')
        diaIdsFile = file(diaIdsFileName, 'r')
        outFile = file(outFileName, 'w')
        
        diaIds = lineNumToIdsAndMjds(diaIdsFile)

        # the track file holds line numbers from the .tracklets.miti
        # file. the diaIds file holds lines s.t. line X of diaIds file
        # is the DiaId and obs. time of the dia source at line X of
        # the .tracklets.miti file.
        
        trackLine = trackFile.readline()
        while trackLine != "":
            lineNums = map(int, trackLine.split())
            trackDiaIds = []
            #trackMjds = []
            for lineNum in lineNums:
                trackDiaIds.append(diaIds[lineNum])
                #trackMjds.append(mjds[lineNum])

            for dia in sorted(trackDiaIds):
                outFile.write(" %d " % dia)

            outFile.write("\n")

            trackLine = trackFile.readline()
            trackId += 1


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

This modified version alters a SINGLE OUTPUT FILE specified as its one argument.

"""



import glob




def lineNumToId(detsFile):
    line = detsFile.readline()
    ids = []
    while line != "":
        diaId = int(line.split()[0])
        ids.append(diaId)        
        line = detsFile.readline()
    return ids











if __name__ == "__main__":

    import sys

    trackFileName = sys.argv[1]
    diaIdsFileName = sys.argv[2]
    outFileName = sys.argv[3]


    print "Reading tracks-by-indices from: ", trackFileName
    print "Reading diaIds from: ", diaIdsFileName
    print "Writing outFile: ", outFileName


    trackFile = file(trackFileName, 'r')
    diaIdsFile = file(diaIdsFileName, 'r')
    outFile = file(outFileName, 'w')
    
    diaIds = lineNumToId(diaIdsFile)

    # the track file holds line numbers from the .tracklets.miti
    # file. the diaIds file holds lines s.t. line X of diaIds file
    # is the DiaId and obs. time of the dia source at line X of
    # the .tracklets.miti file.
    
    trackLine = trackFile.readline()
    while trackLine != "":
        lineNums = map(int, trackLine.split())
        trackDiaIds = []
        for lineNum in lineNums:
            trackDiaIds.append(diaIds[lineNum])
            
        for dia in sorted(trackDiaIds):
            outFile.write(" %d " % dia)
        outFile.write("\n")

        trackLine = trackFile.readline()
    outFile.close()
    print "FINISHED successfully."

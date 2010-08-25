#!/usr/bin/env python

""" jmyers july 13 2010

currently, the MITI format gives only one column for ID; ergo, when
tracklets are created, this ID holds the ID of the tracklet, not the
detection.

As a result, our .tracklet.miti files lose some important information:
the DiaSource IDs of the DiaSources.

This script uses the existing .dets and .ids files and uses this as a
way to determine the DiaSourceIds of each line in the .tracklets.miti files.

"""

import glob


DETS_GLOB=glob.glob("*.dets")
TRACKLETS_SUFFIX=".ids"
OUTFILES_SUFFIX=".tracklets.miti.diaIds"

def lineNumToIdsAndMJDs(detsFile):
    line = detsFile.readline()
    res = []
    while line != "":
        items = line.split()
        diaId = int(items[0])
        mjd = float(items[1])
        res.append([diaId, mjd])
        line = detsFile.readline()
    return res

if __name__=="__main__":

    for detsFile in DETS_GLOB:
        idsFile = detsFile[:-5] + TRACKLETS_SUFFIX
        outFile = detsFile[:-5] + OUTFILES_SUFFIX
        print detsFile, idsFile, outFile

        idsFile = file(idsFile, 'r')
        detsFile = file(detsFile, 'r')
        outFile = file(outFile, 'w')
        
        #get mapping of detection line number to diaId
        detIdsAndMjds = lineNumToIdsAndMJDs(detsFile)
        
        line = idsFile.readline()
        while line != "":
            lineNums = map(int, line.split())
            
            for lineNum in lineNums:
                outFile.write("%d %f\n" % (detIdsAndMjds[lineNum][0], detIdsAndMjds[lineNum][1]) )
                

            line = idsFile.readline()
            

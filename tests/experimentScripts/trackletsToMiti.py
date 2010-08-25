#!/usr/bin/env python

"""

jmyers

may 13 2010

Kubica's version of findTracklets takes input as a set of MITI
detections, where detections with the same ID belong to the same
tracklet.

My version takes a series of MITI detections with unique ids, and a series of indices into
the detections file (NOT DETECTION IDS).

This takes the latter version (created with, say, make20DayWindows.py)
and converts to the former, Kubica version. This will be used for
side-by-side comparison of Kubica C vs. my C++ LinkTracklets versions.

"""


def mitifyTracklets(detsFile, pairsFile, outFile):
    """
    writes tracklets in a linkTracklets-friendly format to outFile    
    using pairsFile and detsFile.  The first ID used is startID; the number 
    of total tracklets written is returned.                                
    """

    trackletsWritten = 0
    allDets = detsFile.readlines()
    for line in pairsFile:
        items = line.split()
        items = map(int, items)
        for item in items:
            det = allDets[item]
            det = det.split()
            #field 0 is the obsID, replace with trackletID.   
            outFile.write("%i " % (trackletsWritten))
            for data in det[1:]:
                outFile.write("%s " % (data))
            outFile.write("\n")
        trackletsWritten += 1


if __name__=="__main__":
    import sys
    detsFile = file(sys.argv[1], 'r')
    idsIn = file(sys.argv[2], 'r')
    mitiOut = file(sys.argv[3], 'w')
    mitifyTracklets(detsFile, idsIn, mitiOut)

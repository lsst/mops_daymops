#!/usr/bin/env python

"""

jmyers oct 13 2010

This script takes a pass over some tracks and finds the unique pairs
of endpoint images - the one from the FIRST diaSource, the one from
the SECOND TO LAST diaSource - and finds the number of total, true,
and false tracks at each unique endpoint image pair.

"""


import sys
import glob

FALSE_DIA_SSM_ID="-1" # the ssmId of a DiaSource which is attributable to non-asteroid sources


def readDias(diasDataFile):

    """ reads a dump of dias, which include diaId expMjd ssmId
    obsHistId for every diaSource.  Returns a dict mapping diaId to
    data."""
    idToDias = {}
    line = diasDataFile.readline()
    while line != "":
        [diaId, expMjd, ssmId, obsHistId] = line.split()
        diaId = int(diaId)
        expMjd = float(expMjd)
        obsHistId = int(obsHistId)
        
        idToDias[diaId] = diaSource(diaId=diaId, obsTime=expMjd, ssmId=ssmId, obsHistId=obsHistId)
        
        line = diasDataFile.readline()

    return idToDias


def lookUpDia(diasLookupDict, diaId):
    return diasLookupDict[diaId]



class diaSource:
    def __init__(self, diaId, obsTime, ssmId, obsHistId):
        self.diaId = diaId
        self.obsTime = obsTime
        self.ssmId = ssmId
        self.obsHistId = obsHistId

    def getDiaId(self):
        return self.diaId

    def getObsTime(self):
        return self.obsTime

    def getSsmId(self):
        return self.ssmId

    def getObsHistId(self):
        return self.obsHistId



def hashIntPair(obsHist1, obsHist2):
    return ",".join(map(str, [obsHist1, obsHist2]))

def intPairFromHash(hash):
    return map(int, hash.split(","))


def countTracksByEndpointImages(diasLookupDict, tracksFile):

    """ pass over the tracks in a file and count the number of
    true/false/total tracks, binned by the image generating each
    ENDPOINT TRACKLET.  I.E. since we only use length-2 tracklets, we
    take the image holding the FIRST diaSource in the track, and the
    image holding the SECOND-TO-LAST diaSource in the track.

    Returns a dictionary mapping obsHistId pair to number of
    total/true/false tracks corresponding to that image pair."""

    trackLine = tracksFile.readline()

    imgPairToCounts = {}
    
    while trackLine != "":
        diaIds = map(int, trackLine.split())
        dias = map(lambda x: lookUpDia(diasLookupDict, x), diaIds)

        ssmIds = set(map(lambda x: x.getSsmId(), dias))

        diasByObsTime = sorted(dias, key=(lambda x: x.getObsTime()))

        obsHist1 = diasByObsTime[0].getObsHistId()
        obsHist2 = diasByObsTime[-2].getObsHistId()
        obsHistHash = hashIntPair(obsHist1, obsHist2)

        if (imgPairToCounts.has_key(obsHistHash)):
            [trueAtImage, falseAtImage] = imgPairToCounts[obsHistHash]
        else:
            [trueAtImage, falseAtImage] = [0, 0]

        #figure out if this is a true track, act accordingly
        if (len(ssmIds) == 1) and (FALSE_DIA_SSM_ID not in ssmIds):
            trueAtImage += 1
        else:
            falseAtImage += 1
                    
        #update imgPairToCounts
        imgPairToCounts[obsHistHash] = [trueAtImage, falseAtImage]
        
        trackLine = tracksFile.readline()

    return imgPairToCounts






def writeImgPairsFile(imgPairToCounts, imgPairCountsOut):
    imgPairCountsOut.write("!obsHistId1 obsHistId2  nTracks_startinghere nTrueTracks_startinghere nFalseTracks_startinghere\n")
    for hash in sorted(obsHistCounts.keys()):
        [nTrue, nFalse] =  imgPairToCounts[hash]
        obsHist1, obsHist2 = intPairFromHash(hash)
        imgPairCountsOut.write("%d %d %d %d %d\n" % (obsHist1, obsHist2, nTrue + nFalse, nTrue, nFalse))


    


if __name__=="__main__":
    import time
    "Starting analysis at ", time.ctime()

    if len(sys.argv) != 3:
        print "USAGE: ", sys.argv[0], " diaDataDump tracksGlobPattern"
        sys.exit(1)
        
    [diaDataDump, tracksGlobPattern] = sys.argv[1:]
    tracksGlob = glob.glob(tracksGlobPattern)

    print "Reading diaSource info from ", diaDataDump
    print "Reading tracks from ", tracksGlob

    diasDataFile = file(diaDataDump,'r')

    print "Reading dump of all Dias at ", time.ctime()
    diasLookupDict = readDias(diasDataFile)
    print "Done. Starting analysis at ", time.ctime()
    t0 = time.time()
    
    for tracks in tracksGlob:

        tracksFile = file(tracks,'r')
        obsHistCounts = countTracksByEndpointImages(diasLookupDict, tracksFile)
        print "Done at ", time.ctime()
        dt = time.time() - t0

        statsOut = tracks + ".stats_byObsHistPair" 
        statsOutFile = file(statsOut,'w')

        print "Writing output to ", statsOut, " starting at ", time.ctime()
        writeImgPairsFile(obsHistCounts, statsOutFile)

        statsOutFile.close()
    
        print "Analysis of ", tracks, " DONE and output written successfully at ", time.ctime()
        
    print " ALL ANALYSES DONE AT ", time.ctime()

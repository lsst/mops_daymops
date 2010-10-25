#!/usr/bin/env python

"""

jmyers oct 18 2010

A script to make a pass over the tracks and count true/false/total
tracks by the cardinality (number of detections linked into the
track) and by the number of nights on which the track was observed.

"""


import sys

FALSE_DIA_SSM_ID="-1" # the ssmId of a DiaSource which is attributable to non-asteroid sources

# I tested this: select max(abs(taiMidPoint - round(taiMidPoint))) from fullerDiaSource;
# returns .4429.
GMT_OFFSET_DAYS=0.0


def mjdToNightNum(mjd):
    #see comments near GMT_OFFSET_DAYS
    return int(mjd)





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





def countTracksByCardinalityAndNumNights(diasLookupDict, tracksFile):

    """ Returns a dictionary mapping cardinality to number of
    total/true/false tracks corresponding in that image pair.  Also
    returns a dictionary mapping number of nights with apparitions to
    a count of total/true/false tracks with apparitions on that many
    nights."""

    trackLine = tracksFile.readline()

    cardinalityToCounts = {}
    numNightsToCounts= {}
    
    while trackLine != "":
        diaIds = map(int, trackLine.split())
        dias = map(lambda x: lookUpDia(diasLookupDict, x), diaIds)

        ssmIds = set(map(lambda x: x.getSsmId(), dias))

        if len(ssmIds) == 1 and (not (FALSE_DIA_SSM_ID in ssmIds)):
            isTrue = True
        else:
            isTrue = False
            
        obsTimes = map(lambda x: x.getObsTime(), dias)
        nightsObserved = set(map(mjdToNightNum, obsTimes))
        numNightsObserved = len(nightsObserved)
        cardinality = len(dias)

        for [dict, key] in [[cardinalityToCounts, cardinality], [numNightsToCounts, numNightsObserved]]:
            if dict.has_key(key):
                [nTrue, nFalse] = dict[key]
            else:
                [nTrue, nFalse] = [0,0]
            if isTrue:
                nTrue += 1
            else:
                nFalse += 1
            dict[key] = [nTrue, nFalse]

        trackLine = tracksFile.readline()

    return cardinalityToCounts, numNightsToCounts






def writeCountsDictToFile(dict, outfile, keyString):
    outfile.write("!%s  nTotalTracks nTrueTracks nFalseTracks\n" % keyString)
    for key in sorted(dict.keys()):
        [nTrue, nFalse] =  dict[key]
        outfile.write("%r %d %d %d\n" % (key, nTrue + nFalse, nTrue, nFalse))


    


if __name__=="__main__":
    import time
    "Starting analysis at ", time.ctime()

    if len(sys.argv) != 5:
        print "USAGE: ", sys.argv[0], " diaDataDump tracksFile cardinalityOutFile numNightsOutfile"
        sys.exit(1)
        
    [diaDataDump, tracks, cardinalityOut, numNightsOut] = sys.argv[1:]

    print "Reading diaSource info from ", diaDataDump
    print "Reading tracks from ", tracks
    print "Pringing stats by image length to ",cardinalityOut
    print "Pringing stats by num nights of observation to ", numNightsOut

    diasDataFile = file(diaDataDump,'r')
    tracksFile = file(tracks,'r')
    cardinalityOutFile = file(cardinalityOut,'w')
    numNightsOutFile = file(numNightsOut,'w')
    

    print "Reading dump of all Dias at ", time.ctime()
    diasLookupDict = readDias(diasDataFile)
    print "Done. Starting analysis at ", time.ctime()
    t0 = time.time()

    cardinalityCounts, numNightsCounts = countTracksByCardinalityAndNumNights(diasLookupDict, tracksFile)
    print "Done at ", time.ctime()
    dt = time.time() - t0
    print "Writing output at ", time.ctime()
    writeCountsDictToFile(cardinalityCounts, cardinalityOutFile, "trackCardinality")
    writeCountsDictToFile(numNightsCounts, numNightsOutFile, "uniqueNights")

    cardinalityOutFile.close()
    numNightsOutFile.close()
    
    print "Analysis DONE and output written successfully at ", time.ctime()

#!/usr/bin/env python

"""

jmyers oct 7 2010

Our mega-analysis script for the output of linkTracklets (check out
that filename!)

The idea is to make a single pass over each .byDiaIds file (tracks,
separated by newlines, each track as expressed as a set of Dia IDs,
separated by whitespace).  Since we tend to pay 'per-pass' in terms of
runtime, this will be a lot faster than running individual scripts to
do each query.

Also needed is a file which maps diaId to ssmId and obsHistId.  This
way we can judge whether a track was true, its corresponding object,
and pair track counts to their start (and possibly someday, end)
images.

So, here are the parameters:

diaDataDump tracksFile statsFile per-obsHistCountsFile foundObjectsFile

the diaDataDump is the file containing (whitespace-delimited fields,
newline-delimited lines):

diaId observingTime(MJD) objectID obsHistId(image where diaId was observed) 



tracksFile will be the tracks from linkTracklets.
-------------------------------------------------------

statsFile will get the following data:

- total tracks in the tracksFile
- true tracks in the tracksFile
- false tracks in the tracksFile
- number of unique findable objects in the tracksFile (objects with some true track)



per-obsHistCountsFile gets the following:
--------------------------------------------------------

for every obsHistId, the number of tracks STARTING IN that image, in the format:

obsHistId numTracks



foundObjectsFile gets the following:
--------------------------------------------------------
the name of every object found (associated with some true track), newline-delimited.

"""


import sys

FALSE_DIA_SSM_ID="-1" # the ssmId of a DiaSource which is attributable to non-asteroid sources

PRELOAD_DIAS_FROM_FILE=False

if not PRELOAD_DIAS_FROM_FILE:
    import MySQLdb as db


    OPSIM_DB="opsim_3_61"
    OPSIM_TABLE="output_opsim3_61"
    
    DIAS_DB="mops_noDeepAstromError"
    DIAS_TABLE="fullerDiaSource"
    
    DB_USER="jmyers"
    DB_PASS="jmyers"
    DB_HOST="localhost"


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



def lookUpDias(diasLookupTool, diaIds):
    if PRELOAD_DIAS_FROM_FILE:
        return map(lambda x: diasLookupTool[x], diaIds)
    else:
        cursor = diasLookupTool
        sql = """ SELECT diaSourceId, taiMidPoint, ssmId, opSimId FROM %s.%s 
                 WHERE diaSourceId IN (""" % (DIAS_DB, DIAS_TABLE)
        first = True;
        for dia in diaIds:
            if not first:
                sql += ", "
            first = False
            sql += str(dia)
        sql += ");"
        cursor.execute(sql)
        rows = cursor.fetchall()
        return map(lambda row: 
                   diaSource(diaId=row[0], obsTime=row[1], ssmId=row[2], obsHistId=row[3]), 
                   rows)


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






def getLotsOfStatsFromTracksFile(diasLookupTool, tracksFile, trueTracksOutFile):

    """return number of true tracks, number of false tracks, a
    dictionary mapping obsHistId (image ID) to true/false track counts
    (as a two-part list), and the set of findable objects

    also, whenever a true track is found, it is written to trueTracksOutFile
    """

    trackLine = tracksFile.readline()

    nTrue = 0
    nFalse = 0
    obsHistCounts = {}
    foundObjects = set()
    lineCount = 0
    
    while trackLine != "":
        lineCount += 1
        print "Read ", lineCount, " lines so far"
        diaIds = map(int, trackLine.split())
        dias = lookUpDias(diasLookupTool, diaIds)

        ssmIds = set(map(lambda x: x.getSsmId(), dias))

        #figure out the first obs time and the obsHistId associated with it
        minObsTime = None
        firstObsHistId = None
        for dia in dias:
            obsTime = dia.getObsTime()
            if minObsTime == None:
                minObsTime = obsTime
                firstObsHistId = dia.getObsHistId()
            elif minObsTime > obsTime:
                minObsTime = obsTime
                firstObsHistId = dia.getObsHistId()

        if obsHistCounts.has_key(firstObsHistId):
            [trueAtImage, falseAtImage] = obsHistCounts[firstObsHistId]
        else:
            trueAtImage = 0
            falseAtImage = 0
            
        #figure out if this is a true track, act accordingly
        if (len(ssmIds) == 1) and (FALSE_DIA_SSM_ID not in ssmIds):
            nTrue += 1
            isTrueTrack = True
            objName = dias[0].ssmId
            foundObjects.add(objName)
            trueAtImage += 1
            trueTracksOutFile.write(trackLine.strip() + '\n')
        else:
            nFalse += 1
            falseAtImage += 1

        #update obsHistCounts
        obsHistCounts[firstObsHistId] = [trueAtImage, falseAtImage]
        
        trackLine = tracksFile.readline()

    return nTrue, nFalse, obsHistCounts, foundObjects




def writeStatsFile(nTrue, nFalse, foundObjects, statsOutFile):
    statsOutFile.write("!num_total_tracks num_true_tracks    num_false_tracks    num_found_objects\n")
    statsOutFile.write("%d %d %d %d\n" % (nTrue+nFalse, nTrue, nFalse, len(foundObjects)))


def writeObsHistFile(obsHistCounts, obsHistCountsOut):
    obsHistCountsOut.write("!obsHistId  nTracks_startinghere nTrueTracks_startinghere nFalseTracks_startinghere\n")
    for obsHistId in sorted(obsHistCounts.keys()):
        [nTrue, nFalse] =  obsHistCounts[obsHistId]
        obsHistCountsOut.write("%d %d %d %d\n" % (obsHistId, nTrue + nFalse, nTrue, nFalse))


def writeFoundObjectsFile(foundObjects, foundObjectsOut):
    for foundObject in foundObjects:
        foundObjectsOut.write("%s\n" % foundObject)
    


if __name__=="__main__":
    import time
    import glob
    "Starting analysis at ", time.ctime()


    # if preloading dias, diasLookupTool is a dict of data read from a
    # file.  if not, diasLookupTool will be a database connection.
    # Reading from a dict is MUCH, MUCH faster (1000x) than making
    # many DB queries.  However, it does impose a few minutes of
    # start-up overhead.

    if PRELOAD_DIAS_FROM_FILE:
        if len(sys.argv) != 3:
            print "USAGE: ", sys.argv[0], " diaDataDump tracksGlob"
            sys.exit(1)
        [diaDataDump, tracksGlobPattern] = sys.argv[1:]
        print "Reading diaSource info from ", diaDataDump
        diasDataFile = file(diaDataDump,'r')
        
        print "Reading dump of all Dias at ", time.ctime()
        diasLookupTool = readDias(diasDataFile)
    else:

        if len(sys.argv) != 2:
            print "USAGE: ", sys.argv[0], " tracksGlob"
            sys.exit(1)

        tracksGlobPattern = sys.argv[1]
        conn = db.connect(user=DB_USER, passwd=DB_PASS, host=DB_HOST)
        diasLookupTool = conn.cursor()

    tracksGlob = glob.glob(tracksGlobPattern)
    

    print "Reading tracks from ", tracksGlob


    for tracks in tracksGlob:
            tracksFile = file(tracks,'r')
            statsOutFile = file(tracks + ".stats",'w')
            obsHistOutFile = file(tracks + ".stats_perObsHist",'w')
            foundObjectsOutFile = file(tracks + ".foundObjects",'w')
            trueTracksOutFile = file(tracks + ".trueTracks.byDiaId", 'w')
            
            print "Starting analysis of ", tracks, " at ", time.ctime()
            t0 = time.time()
            nTrue, nFalse, obsHistCounts, foundObjects = getLotsOfStatsFromTracksFile(diasLookupTool, tracksFile, trueTracksOutFile)
            print "Done at ", time.ctime()
            dt = time.time() - t0
            print "Reading/analyzing ", nTrue + nFalse, " tracks took ", dt, " seconds."
            print "Writing output at ", time.ctime()
            writeStatsFile(nTrue, nFalse, foundObjects, statsOutFile)
            writeObsHistFile(obsHistCounts, obsHistOutFile)
            writeFoundObjectsFile(foundObjects, foundObjectsOutFile)
            
            statsOutFile.close()
            obsHistOutFile.close()
            foundObjectsOutFile.close()
            
            print "Analysis DONE and output written successfully at ", time.ctime()

    print "ALL ANALYSES FINISHED AT ", time.ctime()

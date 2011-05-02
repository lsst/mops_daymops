#!/usr/bin/env python

"""

jmyers oct 14 2010

FUNCTIONALITY:

This script uses the DB to look up the obsHistId of the first
detection in each tracklet, and dump the tracklet to a file with that
obsHistId.

This version makes two passes through the file, one to find all neede
diaSources, then one to do the actual binning.  After the first pass
all needed dias are requested from the DB.  This is in the hopes that
reducing the number of DB accesses may speed things up.

RATIONALE:

In order to get a better load distribution (sorta...) of our
linkTracklets runs, we'd like to distribute workloads by image rather
than night.  To do this we need to separate out the tracklets by first
image time.  Since there are a LOT of tracklets, we're going to put
them into carefully-named files rather than put them into the DB.



USAGE:

<scriptName> trackletsFile outputDirectory
"""

import sys
import os.path
import MySQLdb as db
import mopsDatabases
import time

OUT_TRACKLETS_SUFFIX=".tracklets.byDiaId"


def getAllDiasInFile(inTracklets):
    """ return a set of all diaSources present in tracklets file. """
    line = inTracklets.readline()
    allDias = set()
    while line != "":
        dias = map(int, line.split())
        for d in dias:
            allDias.add(d)
        line = inTracklets.readline()
    return allDias



def getDiaTimesAndImages(dias, dbCurs, diasDb, diasTable):
    diasList = list(dias)
    it = 0
    iterSize = 10000
    curDias = diasList[it*iterSize:(it+1)*iterSize]
    toRet = {}
    while curDias != []:
        s = """ SELECT diaSourceId, taiMidPoint, opSimId 
                FROM  
                  %s.%s""" \
            % (diasDb, diasTable)

        s += """  WHERE diaSourceId IN ( """ 
        first = True
        for d in curDias:
            if not first :
                s+= ", "
            first = False
            s += str(d)
        s += " );"
        dbCurs.execute(s)
        res = dbCurs.fetchall()
        for r in res:
            toRet[r[0]] = [r[1], r[2]]
        it += 1
        curDias = diasList[it*iterSize:(it+1)*iterSize]

    return toRet



def firstObsHistForDias(dias, diasDict):
    """ return the obsHistId of the earliest image represented in the track(let)."""

    if (len(dias)) < 1:
        raise Exception("We need at least one diaSource in arguments.")

    #diasDict maps from diaSourceId -> [obsTime, obsHistId]
    earliestTime, earliestObsHist = diasDict[dias[0]]
    for dia in dias[1:]:
        thisTime, thisObsHist = diasDict[dia]
        if thisTime < earliestTime:
            earliestTime = thisTime
            earliestObsHist = thisObsHist

    return earliestObsHist



def getFileForObsHist(obsHist, outDirectory, obsHistToFile):

    """ looks up the right outfile for the current obsHist, opening it
    if necessary and adding it to obsHistToFile dict."""

    if not obsHistToFile.has_key(obsHist):
        # open the file! save it to dictionary.
        filename = os.path.join(outDirectory, str(obsHist) + OUT_TRACKLETS_SUFFIX)
        newFile = file(filename, 'w')
        obsHistToFile[obsHist] = newFile

    return obsHistToFile[obsHist]
        



def writeTrackletsToPerObsHistFiles(inTracklets, outDirectory, diasDict):
    tletLine = inTracklets.readline()
    obsHistToFile = {}

    while tletLine != "":
        dias = map(int, tletLine.split())
        obsHist = firstObsHistForDias(dias, diasDict)

        outFile = getFileForObsHist(obsHist, outDirectory, obsHistToFile)
        outFile.write("%s\n" % (tletLine.strip()))

        tletLine = inTracklets.readline()
    
    #close all files
    for obsHist in obsHistToFile.keys():
        obsHistToFile[obsHist].close()



if __name__=="__main__":
    print "Reading tracklets from ", sys.argv[1]
    print "Writing output to the directory ", sys.argv[2]
    diasDb = sys.argv[3]
    diasTable = sys.argv[4]
    print "Reading diaSources from DB ", diasDb, ".", diasTable
    inTracklets = file(sys.argv[1],'r')
    outDirectory = sys.argv[2]
    curs = mopsDatabases.getCursor()
    # get all dias in file
    print "Finding all dias in file at ", time.ctime()
    allDias = getAllDiasInFile(inTracklets)
    print "Found ", len(allDias), " dias referenced in file. "
    if (len(allDias)) > 0:
        # fetch those dias from DB.
        print "Fetching needed dias from DB at ", time.ctime()
        diaDict = getDiaTimesAndImages(allDias, curs, diasDb, diasTable)
        print "Read ", len(diaDict.keys()), " dias from DB."
        # return to beginning of file
        inTracklets.seek(0)
        print "Starting to sort tracklets at ", time.ctime()
        writeTrackletsToPerObsHistFiles(inTracklets, outDirectory, diaDict)
    else: 
        print "No processing required due to empty infile."
    print "DONE writing output files at ", time.ctime()
    print ""

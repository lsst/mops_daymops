#!/usr/bin/env python

"""

jmyers March 29 2012

FUNCTIONALITY:

The previous version of this script read tracklets from a single file,
then read diaSources from a DB.

This version reads tracklets from a set of files and fetches
diaSources from a flat file (Lynne format, which includes obsHistId.)
This is advantageous when there are many tracklets and the dias file
is local but the DB is non-local (e.g. running on Teragrid, where the
DB is in Tucson)

It then groups the tracklets by their start image, which greatly
simplifies the process of building inputs for linkTracklets.


USAGE:

<scriptName> trackletsFile diasFile outputDirectory
"""

import sys
import os.path
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



def getDiaTimesAndImages(diasFile):
    toRet = {}
    f = file(diasFile,'r')
    line = f.readline()
    count = 0
    while line != "":
        [diaId, obsHistId, ssmId, ra, decl, mjd, mag, snr] = line.split()
        diaId = int(diaId)
        obsHistId = int(obsHistId)
        mjd = float(mjd)
        toRet[diaId] = [mjd, obsHistId]
        count += 1
        if count % 100000 == 0:
            print "Finished reading ", count, " diaSources so far..."
        line = f.readline()

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
    tFile=sys.argv[1]
    print "Reading tracklets from ", tFile
    diasFile = sys.argv[2]
    print "Reading corresponding Dias (assumed to be in fullerDiaSource format) from ", diasFile
    outDirectory = sys.argv[3]
    print "Writing output to the directory ", outDirectory

    print "Reading diaSources from file ", diasFile, " at ", time.asctime()
    diaDict = getDiaTimesAndImages(diasFile)
    print "Finished reading diaSources from file ", diasFile, " at ", time.asctime()

    print "Reading tracklets in ", tFile, " and sorting, starting at ", time.asctime()
    inTracklets = file(tFile,'r')
    writeTrackletsToPerObsHistFiles(inTracklets, outDirectory, diaDict)
    inTracklets.close()
    print "Finished processing tracklets in ", tFile, " at ", time.asctime()
    print ""

#!/usr/bin/env python

"""

jmyers oct 14 2010

FUNCTIONALITY:

This script uses the DB to look up the obsHistId of the first
detection in each tracklet, and dump the tracklet to a file with that
obsHistId.


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
import time

OPSIM_DB="opsim_3_61"
OPSIM_TABLE="output_opsim3_61"

DIAS_DB="mops_noDeepAstromError"
DIAS_TABLE="fullerDiaSource"

DB_USER="jmyers"
DB_PASS="jmyers"

OUT_TRACKLETS_SUFFIX=".tracklets.byDiaId"



def firstObsHistForDias(dias, dbCurs):
    s = """ SELECT ops.obsHistId, ops.expMjd 
            FROM  
              %s.%s AS ops
            JOIN
              %s.%s as dias 
            ON 
             dias.opSimId=ops.obsHistId""" \
        % (OPSIM_DB, OPSIM_TABLE, DIAS_DB, DIAS_TABLE)
        
    s += """  WHERE dias.diaSourceId IN ( """ 
    for i in range(len(dias)):
        s += str(dias[i])
        if i < len(dias) -1 :
            s+= ","

    s += """)
            GROUP BY ops.expMjd 
            ORDER BY ops.expMjd
            ;"""
    dbCurs.execute(s)
    res = dbCurs.fetchall()
    # we get the earliest (by date) image first, and the first column is obsHistId.
    return res[0][0]



def getFileForObsHist(obsHist, outDirectory, obsHistToFile):

    """ looks up the right outfile for the current obsHist, opening it
    if necessary and adding it to obsHistToFile dict."""

    if not obsHistToFile.has_key(obsHist):
        # open the file! save it to dictionary.
        filename = os.path.join(outDirectory, str(obsHist) + OUT_TRACKLETS_SUFFIX)
        newFile = file(filename, 'w')
        obsHistToFile[obsHist] = newFile

    return obsHistToFile[obsHist]
        



def writeTrackletsToPerObsHistFiles(inTracklets, outDirectory, dbCurs):
    tletLine = inTracklets.readline()
    obsHistToFile = {}

    while tletLine != "":
        dias = map(int, tletLine.split())
        obsHist = firstObsHistForDias(dias, dbCurs)

        outFile = getFileForObsHist(obsHist, outDirectory, obsHistToFile)
        outFile.write("%s\n" % (tletLine.strip()))

        tletLine = inTracklets.readline()
    
    #close all files
    for obsHist in obsHistToFile.keys():
        obsHistToFile[obsHist].close()



if __name__=="__main__":
    inTracklets = file(sys.argv[1],'r')
    outDirectory = sys.argv[2]
    conn = db.connect(user=DB_USER, passwd=DB_PASS)
    curs = conn.cursor()
    print "Starting to sort tracklets at ", time.ctime()
    writeTrackletsToPerObsHistFiles(inTracklets, outDirectory, curs)
    print "DONE writing output files at ", time.ctime()

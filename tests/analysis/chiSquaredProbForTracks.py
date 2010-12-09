#!/usr/bin/env python
"""
Given a .c.tracks.byDiaIds file, get the chiSquaredProbability of tracks within it

Adapted from Tim Axelrod's AnalyzeTracks6.py

"""
import string, sys, os
import glob
import re
import math
import numpy as np
import MySQLdb
#import matplotlib.pyplot as plt
import numpy as np
import scipy.stats.stats as st

# set up mysql access params
mySqlHost = 'localhost'
mySqlUser = 'jmyers'
mySqlPasswd = 'jmyers'
mySqlDb = 'mops_noDeepAstromError'
mySqlEphemDb='dc3b_ephem_nodeep'
    
db = MySQLdb.connect(host=mySqlHost, user=mySqlUser, passwd=mySqlPasswd, db=mySqlDb)
ephemDb = MySQLdb.connect(host=mySqlHost, user=mySqlUser, passwd=mySqlPasswd, db=mySqlEphemDb)
dtMin = 0.01
nominalAstroErr = 0.1/3600.0





def getSsmTrack(ssmId):
    trackDia=[]
    trackT=[]
    trackRa=[]
    trackDec=[]

    c=ephemDb.cursor()
    query = 'SELECT diaSourceID, taiMidPoint, ra, decl from dias_nodeep where ssmID=%s order by taiMidPoint' % (ssmId)
    c.execute(query)
    result = c.fetchall()
    for row in result:
        trackDia.append(row[0])
        trackT.append(row[1])
        trackRa.append(row[2])
        trackDec.append(row[3])

    return (trackDia,trackT,trackRa,trackDec)
        
def getTrackDias(diaList):
    trackT=[]
    trackRa=[]
    trackDec=[]
    trackSsmId=[]
    c=db.cursor()

    for dia in diaList:
        query = "select taiMidPoint,ra,decl,ssmId from fullerDiaSource where diaSourceId=%d" % (int(dia))
        c.execute(query)
    
        result = c.fetchall()
        trackT.append(result[0][0])
        trackRa.append(result[0][1])
        trackDec.append(result[0][2])
        trackSsmId.append(result[0][3])

    c.close()
    return (trackT,trackRa,trackDec,trackSsmId)

def analyzeTrack(diaList):

    nTrackPts = len(diaList)

    (t, ra, dec, ssmId)= getTrackDias(diaList)

# check to see whether this is a real track

    ssmId = np.array(ssmId)
    bad = np.where(ssmId != ssmId[0])
    if len(bad[0]):
        trackOK = False
    else:
        trackOK = True
# fit quadratics in t to ra and dec

    raFit = np.polyfit(t,ra, 2)
    raFitPts = np.polyval(raFit, t)
    raChisq = np.sum((raFitPts - ra)**2)/nominalAstroErr**2
    raChisqProb = st.chisqprob(raChisq, nTrackPts)
    
    decFit = np.polyfit(t,dec, 2)
    decFitPts = np.polyval(decFit, t)
    decChisq = np.sum((decFitPts - dec)**2)/nominalAstroErr**2
    decChisqProb = st.chisqprob(decChisq, nTrackPts)
    
#    print raChisq, decChisq, raChisqProb, decChisqProb
    
# calculate the chisq for each fit, assuming for the moment, fixed error
# in each ra and dec measurement of 0.1"

# calculate the chisq probability

    
    return (trackOK, ssmId[0], t, ra, dec, raChisqProb, decChisqProb)


def analyzeTrackSet(fileNameIn, fileNameOut, i1=None, i2=None, plot=False):

    trackFile = open(fileNameIn, 'r')
    outFile = open(fileNameOut, 'w')
    print trackFile, outFile
    tracks = trackFile.readlines()

    print "Read ", len(tracks), " tracks into memory."

    outFile.write("!track_is_true ssmId ra+DecChiSqProb raChiSqProb decChiSqProb [track dias]\n")

    if plot:
        plt.hold(1)

    if i1 == None:
        i1 = 0
    if i2 == None:
        i2 = len(tracks)

    outCount = 0
    for trackIndex in range(i1, i2):
        if outCount % 1000 == 0:
            print "Working on track ", outCount, " of ", i2 - i1
        track = tracks[trackIndex]
        trackDias = track.split()
        (trackOK, ssmId, t, ra, dec, raChisqProb, decChisqProb) = analyzeTrack(trackDias)
        if trackOK:
            ssmId = str(ssmId)
        else:
            ssmId = "FalseTrack"
        print >>outFile, '%r %s %e %e %e [%s]' % \
            (trackOK, ssmId, raChisqProb+decChisqProb, raChisqProb, decChisqProb, track.strip())
        if plot:
            plt.figure()
            plotTrack( t, ra, dec)
            plt.draw()
            response = raw_input('continue?')
            if response=='n':
                outFile.close()
                return
    
        outCount += 1

if __name__=="__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]

    analyzeTrackSet(infile, outfile)

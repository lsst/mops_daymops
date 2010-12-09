#!/usr/bin/python

"""

jmyers

Given a set of tracklets as diaIds, find the one true tracklet per
object seen which has the longest deltaTime from first to last
endpoint, then write this to the output file.  False tracklets are all
written to output.


The original impetus for this file was to be used in order to
understand the relationship between redundant true tracklets and track
output set size.  By removing all but the best tracklet associated
with each object, we expect the number of tracks to decrease
dramatically.

"""


import MySQLdb
import string, sys, os

# set up mysql access params
mySqlHost = 'localhost'
mySqlUser = 'jmyers'
mySqlPasswd = 'jmyers'
mySqlDb = 'mops_noDeepAstromError'
    
db = MySQLdb.connect(host=mySqlHost, user=mySqlUser, passwd=mySqlPasswd, db=mySqlDb)



def getTrackletDias(diaIds, dbConn=db):
    trackDia=[]
    trackT=[]
    trackSsmId=[]

    c=dbConn.cursor()
    query = 'SELECT diaSourceID, taiMidPoint, ssmId from fullerDiaSource where diaSourceId IN (' 
    for i in range(len(diaIds)):
        query += str(diaIds[i])
        if i < len(diaIds) - 1:
            query += ", "
    
    query += ') order by taiMidPoint'
    c.execute(query)
    result = c.fetchall()
    for row in result:
        trackDia.append(row[0])
        trackT.append(row[1])
        trackSsmId.append(row[2])

    return (trackDia,trackT,trackSsmId)



def getBestTrackletPerObj(tracklets):

    # should be a map from ssmId to [[trackletDias], trackletDt]
    bestTrackletPerObj = {}
    bestTracklets = []
    for tracklet in tracklets:

        (dias, times, ssmIds) = getTrackletDias(tracklet)

        if len(set(ssmIds)) == 1 and (not -1 in ssmIds) and (not "NS" in ssmIds):
            # this is a true tracklet
            ssmId = ssmIds[0]
            trackletDt = times[-1] - times[0]
            bestTrackletDt = -1.
            if bestTrackletPerObj.has_key(ssmId):
                bestTrackletDt = bestTrackletPerObj[ssmId][1]                
            if trackletDt > bestTrackletDt:
                bestTrackletPerObj[ssmId] = [tracklet, trackletDt]
        else:

            # false tracklets pass through this filter happily.
            bestTracklets.append(tracklet)

    # collect up the best true tracklets per object
    for obj in bestTrackletPerObj.keys():
        bestTracklets.append(bestTrackletPerObj[obj][0])

    return bestTracklets


if __name__=="__main__":
    infile = file(sys.argv[1], 'r')
    print "Reading input from ", infile
    tracklets = infile.readlines()
    tracklets = map(lambda x: map(int, x.split()), tracklets)
    print "Done. opening output file at ", sys.argv[2]
    outfile = file(sys.argv[2],'w')
    print "Finding best tracklet per object seen."
    bestTracklets = getBestTrackletPerObj(tracklets)
    print "Done. Writing output."
    for tracklet in bestTracklets:
        for diaId in tracklet:
            outfile.write("%d " % diaId)
        outfile.write("\n")
    outfile.close()
    print "Done writing output."

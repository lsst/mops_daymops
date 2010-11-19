#!/usr/bin/env python


""" jmyers nov 19

For each detection d1 in all tracklets:

   Find the number of unique images holding detections with which d1 is linked.

This is intended to help determine whether we are indeed getting
"stealth" deep stacks due to image overlaps and will likely not be
useful for anything else.

Takes the name of a tracklets file. Uses database to find second endpoints.
"""

import MySQLdb as db
import sys
import glob


DB_USER="jmyers"
DB_PASSWD="jmyers"
DB_HOST="localhost"

OPSIM_DB="opsim_3_61"
OPSIM_TABLE="output_opsim3_61"

DIAS_DB="mops_noDeepAstromError"
DIAS_TABLE="fullerDiaSource"


def getObsHist(diaId, cursor):
    s = """ SELECT dias.opSimId FROM %s.%s as dias """ % (DIAS_DB, DIAS_TABLE)
    s += " WHERE dias.diaSourceId = %d;"""% (diaId)
    cursor.execute(s)
    return cursor.fetchall()[0][0]



def countTrackletSecondEndpointImagesPerDetection(trackletFile, cursor):
    line = trackletFile.readline()
    detIdToImages = {}
    while line != "":
        allDias = map(int, line.split())
        firstDia = allDias[0]
        image = getObsHist(allDias[-1], cursor)
        if (detIdToImages.has_key(firstDia)):
            detIdToImages[firstDia].add(image)
        else:
            detIdToImages[firstDia] = set([image])
            
        line = trackletFile.readline()
    detIdToNumImages = {}
    for detId in detIdToImages:
        #print "detection ", detId, " linked with images ", detIdToImages[detId]
        detIdToNumImages[detId] = len(detIdToImages[detId])

    return detIdToNumImages


if __name__=="__main__":

    trackletFile = file(sys.argv[1],'r')
    conn = db.connect(host=DB_HOST, user=DB_USER, passwd=DB_PASSWD)
    curs = conn.cursor()

    histDict = countTrackletSecondEndpointImagesPerDetection(trackletFile, curs)
    
    print "! diaId numSecondEndpointImages"
    for detId in histDict:
        print detId, histDict[detId]

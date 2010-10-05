#!/usr/bin/env/python

""" jmyers

Calculate RA and Dec initial and changes for tracks we find.  Also
while we're at it check whether the first two and last two detections
were correctly linked.

""" 

import numpy 
import sys
import MySQLdb as db

DIAS_DB="mops_noDeepAstromError"
DIAS_TABLE="fullerDiaSource"

DB_HOST="localhost"
DB_USER="jmyers"
DB_PASS="jmyers"

# set to something positive and we will only read the first
# MAX_TO_READ items from the file.
MAX_TO_READ=-1



figformat = "PS"


def getRaDecObjNameFromDias(dias, cursor):
    sql = """ SELECT ra, decl, ssmId FROM %s.%s WHERE diaSourceId in (""" %\
        (DIAS_DB, DIAS_TABLE)
    firstOne = True

    for dia in dias:
        if not firstOne:
            sql += ", "
        else:
            firstOne = False
        sql += str(dia)

    sql += ");"
    cursor.execute(sql)
    return cursor.fetchall()

def angularDistance(a, b):

    """ return distance between a and b, where a and b are angles in
    degrees. """

    while abs(a - b) > 180:
        if a > b:
            b += 360.
        else:
            a += 360.
    return a - b


def makeContiguous(angles):
    """ given a set of angles (say, RAs or Decs of observation) which
    span a fairly short arc but may actually cross the 0/360 line,
    make these contiguous by using negative angles or whatever is
    necessary.  if this set of angles does NOT span a short arc (>180
    deg) expect all hell to break loose."""
    a0 = angles[0]
    output = [a0]
    for angle in angles[1:]:
        while abs(angle - a0) > 180:
            if angle > a0:
                angle -= 360.
            else:
                angle += 360.
        output.append(angle)
    return output


if __name__=="__main__":
    tracksFileName = sys.argv[1]    
    print "Reading ", tracksFileName, " as a set of sets of DiaIds..."
    tracksFile = file(tracksFileName,'r')
    trueTrackStatsName = sys.argv[2]
    falseTrackStatsName = sys.argv[3]
    print "Writing true track stats to ", trueTrackStatsName
    print "Writing false track stats to", falseTrackStatsName
    

    outTrue = file(trueTrackStatsName, 'w')
    outFalse = file(falseTrackStatsName,'w')

    header="""!ra0 dec0 dRa dDec trackTrue firstTrackletTrue lastTrackletTrue\n"""
    outTrue.write(header)
    outFalse.write(header)

    dbConn = db.connect(db=DIAS_DB, user=DB_USER, host=DB_HOST, passwd=DB_PASS)
    dbCurs = dbConn.cursor()

    lineNum = 0
    line = tracksFile.readline()

    while line != "":
        
        dias = map(int, line.split())
        
        rasDecsNames = getRaDecObjNameFromDias(dias, dbCurs)

        ras   = map(lambda x: x[0], rasDecsNames)
        decs  = map(lambda x: x[1], rasDecsNames)
        names = map(lambda x: x[2], rasDecsNames)
        
        ras = makeContiguous(ras)
        decs = makeContiguous(decs)
        dRa = angularDistance(ras[0], ras[-1])
        dDec = angularDistance(decs[0], decs[-1])

        if len(set(names)) == 1:
            isTrue = True
        else:
            isTrue = False

        #this doesn't really deal with length > 2 tracklets but close enough
        firstTrackletTrue = (len(set(names[:2])) == 1)
        lastTrackletTrue = (len(set(names[-2:])) == 1)

        #print ras
        #print decs
        #print names
        #print dRa, dDec, isTrue, firstTrackletTrue, lastTrackletTrue

        outStr = "%10f %10f %10f %10f %r %r %r\n" % \
            ( ras[0], decs[0], dRa, dDec, isTrue, 
              firstTrackletTrue, lastTrackletTrue)

        if isTrue:
            outTrue.write(outStr)
        else:
            outFalse.write(outStr)


        if (MAX_TO_READ > 0) and (MAX_TO_READ < lineNum):
            print "Read ", lineNum, " items, stopping"
            line = ""
        else:
            line = tracksFile.readline()
            lineNum += 1

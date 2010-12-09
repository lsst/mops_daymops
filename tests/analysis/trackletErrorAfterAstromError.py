#!/usr/bin/env python

""" jmyers dec 1 2010

We're currently interested in the velocity errors on tracklets caused
by astrometric error applied to detections (relative to the "true"
velocity given by the ephemerides without astrometric error).
Particularly, we are interested in how this relates to the time span
of the tracklet.

This script should take a tracklets file (holding tracklets as sets of
diaSource Ids) and will calculate the ra/dec velocity using the
first/last detections in the tracklet, using both the pre- and post-
astrometric error diaSources from a database.

Output should be:

deltaTime firstDetectionRaError firstDetectionDecError lastDetectionRaError lastDetectionRaError raVelocityError decVelocityError

one per line, one per tracklet, in the same order as the input tracklets.

"""



import sys
import MySQLdb as db

PUREEPHEMDB="dc3b_ephem_nodeep"
PUREEPHEMTABLE="dias_nodeep"

ASTROMERRDB="mops_noDeepAstromError"
ASTROMERRTABLE="fullerDiaSource"

DB_USER="jmyers"
DB_PASS="jmyers"
DB_HOST="localhost"



def getDiaRaDecMjd(diaId, table, cursor):
    s = """ SELECT ra, decl, taiMidPoint FROM %s WHERE diaSourceId=%d; """ % \
        (table, diaId)
    cursor.execute(s)

    # should raise an exception if the results aren't formatted as
    # expected, or we get multiple hits.  This is a good thing.
    [[ra, dec, mjd]] = cursor.fetchall()
    return [ra,dec, mjd]


def getDtAndErrorsForTracklet(firstDia, lastDia, cursor):
    pureEphemTable = PUREEPHEMDB + "." + PUREEPHEMTABLE
    
    noErrFirst  = getDiaRaDecMjd(firstDia, pureEphemTable, cursor)
    noErrLast   = getDiaRaDecMjd(lastDia,  pureEphemTable, cursor)
    
    withAstromErrTable = ASTROMERRDB + "." + ASTROMERRTABLE
    
    withErrFirst = getDiaRaDecMjd(firstDia, withAstromErrTable, cursor)
    withErrLast  = getDiaRaDecMjd(lastDia,  withAstromErrTable, cursor)

    dt = noErrLast[2] - noErrFirst[2]

    firstDetErrRa = withErrFirst[0] - noErrFirst[0]
    firstDetErrDec = withErrFirst[1] - noErrFirst[1]
    lastDetErrRa = withErrLast[0] - noErrLast[0]
    lastDetErrDec = withErrLast[1] - noErrLast[1]

    noErrRaV =  (noErrLast[0] - noErrFirst[0]) / dt
    noErrDecV = (noErrLast[1] - noErrFirst[1]) / dt
    withErrRaV =  (withErrLast[0] - withErrFirst[0]) / dt
    withErrDecV = (withErrLast[1] - withErrFirst[1]) / dt
    
    return [dt, 
            [firstDetErrRa, firstDetErrDec],
            [lastDetErrRa, lastDetErrDec],
            [withErrRaV - noErrRaV, withErrDecV - noErrDecV]]
    

if __name__=="__main__":
    trackletsIn = file(sys.argv[1],'r')
    outFile = file(sys.argv[2],'w')
    conn = db.connect(host=DB_HOST, user=DB_USER, passwd=DB_PASS)
    curs = conn.cursor()
    
    outFile.write("! trackletDt firstDetRaErr firstDetDecErr lastDetRaErr lastDetDecErr raVelocityErr decVelocityErr\n")

    line = trackletsIn.readline()
    while line != "":
        dias = map(int, line.split())
        [dt, firstDetErrs, lastDetErrs, velocityErrs] = getDtAndErrorsForTracklet(dias[0], dias[-1], curs)

        outFile.write("%10e %10e %10e %10e %10e %10e %10e\n" % (dt, firstDetErrs[0], firstDetErrs[1], 
                                                                lastDetErrs[0], lastDetErrs[1], 
                                                                velocityErrs[0], velocityErrs[1]))
        line = trackletsIn.readline()

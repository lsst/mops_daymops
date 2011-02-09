#!/usr/bin/env python

""" jmyers feb 8 2011


Takes a massive set of tracks (byDiaIds) and uses the database of all
diaSources to find IDS of diaSources NOT present in any true track.

Takes as argument a series of tracks files, or a glob, e.g.

findUnattributedDetections.py run1.byDiaIds run2.byDiaIds

or

findUnattributedDetections.py \*.byDiaIds


Results will be written to stdout.
"""


DB_USER="jmyers"
DB_PASS="jmyers"
DB_DB="mops_noDeepAstromError"
DB_TABLE="fullerDiaSource"

# secret key for a non-asteroid detection.  any track containing this
# a diaSource with this SSM ID is assumed to be false.
FALSE_DET_SSMID=-1


import MySQLdb as db
import glob, sys



def getAllDiaIds(curs):
    query = """ SELECT DISTINCT(diaSourceId) FROM %s.%s; """ % \
        (DB_DB, DB_TABLE)
    curs.execute(query)
    rows = curs.fetchall()
    allDias = map(lambda x: x[0], rows)
    return allDias


def isTrueTrack(curs, dias):
    query = """ SELECT DISTINCT(ssmId) FROM %s.%s 
                   WHERE diaSourceId IN (""" % \
        (DB_DB, DB_TABLE)
    first = True
    for dia in dias:
        if not first:
            query += ", "
        query += " %d " % dia
        first = False
    query += ");"
    curs.execute(query)
    ssmIds = map(lambda x: x[0], curs.fetchall())
    if (len(set(ssmIds)) == 1) and (FALSE_DET_SSMID not in ssmIds):
        return True
    else:
        return False



def findUnattributedDetections(curs, inFiles):
    """ curs should be a DB cursor. inFiles should be a list of open
    files.  curs.DB_TABLE will be searched for all DiaSources, then
    inFiles will be examined and true tracks found. DiaSources NOT in
    some true track are returned."""
    
    #find all dia IDs.
    sys.stderr.write("Reading all diaSource IDs from DB.\n")
    allDias = set(getAllDiaIds(curs))
    sys.stderr.write("Done reading diaSourceIDs from DB.\n")

    for inf in inFiles:
        sys.stderr.write("reading contents of infile %r\n" % inf)
        track = inf.readline()
        while track != "":
            dias = map(int, track.split())
            if isTrueTrack(curs, dias):
                allDias = allDias.difference(set(dias))
            track = inf.readline()

    return allDias



if __name__=="__main__":
    db = db.connect(user=DB_USER, passwd=DB_PASS, db=DB_DB)
    curs = db.cursor()
    
    inFiles = []
    # read each argument. It may be a glob. 
    for arg in sys.argv[1:]: 
        newInFiles = glob.glob(arg)
        # open each file to make sure it exists.
        for inf in newInFiles:
            sys.stderr.write("Opening %r as input file\n" % inf)
            inFiles.append(file(inf,'r'))
            sys.stderr.write("Opened it \n")

    #sys.stderr.write("Got input files: %r\n" % inFiles)
    sys.stderr.write("Calling findUnattributedDetections\n")
    unAttributed = findUnattributedDetections(curs, inFiles)
    sys.stderr.write("Writing results to stdout.")
    for a in unattributed:
        sys.stdout.write("%i\n" % a)
    sys.stdout.flush()
    sys.stderr.write("Done successfully.")

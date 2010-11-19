#!/usr/bin/env python

""" jmyers nov 19 2010

Takes tracklets (or tracks) in the by-dia-id format, looks them up in
the DB, and then counts the number of field crossers and
non-field-crossers. Reports to stdout.

"""

import os.path
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



def countUniqueFields(dias, cursor):
    s = """ SELECT COUNT(DISTINCT(ops.fieldId)) FROM %s.%s AS ops
                JOIN %s.%s as dias
            ON
               dias.opsimId=ops.obsHistId """ % (OPSIM_DB, OPSIM_TABLE, DIAS_DB, DIAS_TABLE)
    s += " WHERE dias.diaSourceId IN ( """
    for i in range(len(dias)):
        if (i != 0):
            s += ", "
        s += str(dias[i])
    s += ");"
    cursor.execute(s)
    return cursor.fetchall()[0][0]

def countFieldCrossers(trackletsFile, cursor):
    """ takes an open file of tracklets (or tracks) and a DB
    cursor. Counts the number of field crossers (and number of
    non-field crossers) and returns both in that order."""

    line = trackletsFile.readline()
    fieldCrossers = 0
    sameFielders = 0
    while line != "":
        dias = map(int, line.split())
        nFields = countUniqueFields(dias, cursor)
        if nFields == 1:
            sameFielders += 1
        else:
            fieldCrossers += 1
        line = trackletsFile.readline()
    return fieldCrossers, sameFielders





if __name__=="__main__":
    trackletsFile = file(sys.argv[1], 'r')
    conn = db.connect(host=DB_HOST, user=DB_USER, passwd=DB_PASSWD)
    curs = conn.cursor()
    
    fieldCrossers, sameFielders = countFieldCrossers(trackletsFile, curs)
    
    print "! Infile nFieldCrossers nSameFielders nTotal"
    print sys.argv[1], fieldCrossers, sameFielders, fieldCrossers + sameFielders

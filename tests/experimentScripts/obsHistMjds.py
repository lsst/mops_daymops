#!/usr/bin/env python

""" jmyers oct 8

For every one of a variety of obsHistIds, count how many times that field was revisited within 15 days.

"""

import sys
import MySQLdb as db


#in days. We return only revisits which occur withing TIME_WINDOW of the requested image.
TIME_WINDOW=15.0

DB_USER="jmyers"
DB_PASS="jmyers"
DB_HOST="localhost"

OPSIM_DB="opsim_3_61"
OPSIM_TABLE="output_opsim3_61"





def timeOfObsHist(obsHistId, cursor):
    s = """ SELECT expMJD FROM %s.%s 
            WHERE obsHistId=%d
            ;""" \
    % (OPSIM_DB, OPSIM_TABLE, obsHistId)
    #print s
    cursor.execute(s)
    res = cursor.fetchall()[0][0]
    return res

if __name__=="__main__":
    inObsHistsFile=file(sys.argv[1],'r')

    obsHists = map(int, inObsHistsFile.readlines())

    conn = db.connect(host=DB_HOST, user=DB_USER, passwd=DB_PASS)
    cursor = conn.cursor()

    for i in obsHists:
        first = timeOfObsHist(i, cursor)
        print i, first

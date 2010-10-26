#!/usr/bin/env python


""" jmyers oct 8


Reads through opSim and returns the number of Images which contained
DiaSources and happened between some min and max time.

"""

DB_USER="jmyers"
DB_PASS="jmyers"
DB_HOST="localhost"

DIAS_DB="mops_noDeepAstromError"
DIAS_TABLE="fullerDiaSource"
OPSIM_DB="opsim_3_61"
OPSIM_TABLE="output_opsim3_61"

import MySQLdb as db
import sys


def listImagesBetween(cursor, dateA, dateB, epsilon=1e-5):
    s = """select ops2.obsHistId, ops2.expMjd, ops2.fieldId 
            from
          %s.%s as ops2

             join (
                 select distinct(ops.obsHistId) as obsHistId
                 from %s.%s as ops 
                  join 
                 %s.%s as dia 
                   on ops.obsHistId = dia.opSimId
                 where ops.expMJD > %f and ops.expMjd < %f
                 group by ops.expMjd              
                 order by ops.expMjd) as subq

            on subq.obsHistId = ops2.obsHistId
            ;""" \
        % (OPSIM_DB, OPSIM_TABLE, OPSIM_DB, OPSIM_TABLE, DIAS_DB, DIAS_TABLE, dateA - epsilon, dateB + epsilon)
    cursor.execute(s)
    res = cursor.fetchall()
    return res




if __name__=="__main__":
    firstDate = float(sys.argv[1])
    lastDate = float(sys.argv[2])

    conn = db.connect(host=DB_HOST, user=DB_USER, passwd=DB_PASS)
    curs = conn.cursor()

    imgs = listImagesBetween(curs, firstDate, lastDate) 
    for img in imgs:
        print img[0], img[1], img[2]
                        
    

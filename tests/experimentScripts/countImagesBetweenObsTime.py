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


def numImagesBetween(cursor, dateA, dateB, epsilon=1e-5):
    s = """select count(distinct(ops.expMJD)) 
            from %s.%s as ops 
              join 
            %s.%s as dia 
             on ops.obsHistId = dia.opSimId
            where ops.expMJD > %f and ops.expMjd < %f
            ;""" \
        % (OPSIM_DB, OPSIM_TABLE, DIAS_DB, DIAS_TABLE, dateA - epsilon, dateB + epsilon)
    cursor.execute(s)
    res = cursor.fetchall()
    if (len(res)) != 1:
        raise Exception("only expected one results row, got %d" % (len(res)))
    if (len(res[0])) != 1:
        raise Exception("only expected one results col, got %d" % (len(res)))

    return res[0][0]

def numDiasBetween(cursor, dateA, dateB, epsilon=1e-5):
    s = """select count(*) 
            from 
            %s.%s as dia 
            where dia.taiMidPoint > %f and dia.taiMidPoint < %f
            ;""" \
        % (DIAS_DB, DIAS_TABLE, dateA - epsilon, dateB + epsilon)
    cursor.execute(s)
    res = cursor.fetchall()
    if (len(res)) != 1:
        raise Exception("only expected one results row, got %d" % (len(res)))
    if (len(res[0])) != 1:
        raise Exception("only expected one results col, got %d" % (len(res)))

    return res[0][0]



if __name__=="__main__":
    firstDate = float(sys.argv[1])
    lastDate = float(sys.argv[2])

    conn = db.connect(host=DB_HOST, user=DB_USER, passwd=DB_PASS)
    curs = conn.cursor()

    imgs = numImagesBetween(curs, firstDate, lastDate) 
    dias = numDiasBetween(curs, firstDate, lastDate)
    print firstDate, lastDate, imgs, dias
                        
    

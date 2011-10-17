#!/usr/bin/env python

import MySQLdb as db
import sys
from glob import glob



def expMjdForFile(lines):
    for line in lines:
        if line[0] == "#":
            items = line.split()
            if len(items)==3:
                if items[1] == "Opsim_expmjd":
                    return float(items[2])
    raise Exception("Could not find line like \"# Opsim_expmjd <num>\"?!")



def addToDb(cursor, diaId, opsimId, ssmId, ra, decl, expMjd, mag):
    sql = """ INSERT INTO fullSkyOneMonth.mopsDetections VALUES
(%d, %d, %d, %f, %f, %f, %f, NULL); """ % \
        (diaId, opsimId, ssmId, ra, decl, expMjd, mag)
    cursor.execute(sql)
    

def getImageLimitingMag(cursor, opsimId):
    sql = """
SELECT 5sigma_modified FROM opsim_3_61.output_opsim3_61 
   WHERE obsHistId = %d;""" % opsimId
    cursor.execute(sql)
    res = cursor.fetchall()[0][0]
    return res


if __name__=="__main__":
    conn = db.connect(user="jmyers", passwd="jmyers", host="localhost")
    cursor = conn.cursor()
    
    #cursor.execute("CREATE DATABASE fullSkyOneMonth;");
    #cursor.execute("CREATE TABLE fullSkyOneMonth.mopsDetections LIKE mops_noDeepAstromError.fullerDiaSource;")

    trimFiles = sorted(glob("*/*.dat"))
    chkpt = file("finishedIngests.txt",'w')
    curDiaId = 0
    
    aboveLimits = file("tooFaint.dets",'w')
    aboveLimits.write("#dummy opSimId ssmId ra decl expMjd mag NULL\n")

    for trimFile in trimFiles:
        # the following line works for this data dump but be careful if reusing this code...
        opsimId = int(trimFile[-12:-4])
        imageLimit = getImageLimitingMag(cursor, opsimId)
        
        lines = file(trimFile,'r').readlines()
        expMjd = expMjdForFile(lines)
        for line in lines:
            if line[0] != "#":
                items = line.split()
                ssmId = int(items[0])
                ra,dec  = map(float, items[1:3])
                mag = float(items[5])
                 
                
                if mag <= imageLimit:
                    addToDb(cursor, 
                            curDiaId, opsimId, ssmId, 
                            ra, dec, expMjd, mag)
                    
                    curDiaId += 1
                else:
                    aboveLimits.write("0 %d %d %d %f %f %f %f NULL\n" % \
                                      (curDiaId, opsimId, ssmId, 
                                       ra, dec, expMjd, mag))

        # checkpoint that we finished this file
        print "Finished ingesting ", trimFile
        chkpt.write("%s %d\n" % (trimFile, curDiaId))


    # done!

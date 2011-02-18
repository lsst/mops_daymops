#!/usr/bin/env python

""" jmyers feb 18 2011

convert tracks (or tracklets even) from byDiaId format (one line per
track, white-space delimited DiaIds which build the track) into
something like the format of Lynne's original DiaSource dump:

#!trackNum diaSourceId opSimId ssmId ra decl taiMidPoint mag snr

Tracks will be assigned IDs sequentially starting at zero.

Since the other version was too slow, this version reads the whole
diaSource DB into memory at start-up.
"""



import mopsDatabases
import sys
import MySQLdb as db





def getDiaInfo(diaIds, idToDataMap):
    return map(lambda x: idToDataMap[x], diaIds)


def fetchAllDiasFromDb(cursor):
    """ return a dictionary mapping diaId to diaSource, for all Dias
    in our data set."""
    toRet = {}
    print "Fetching DiaSources from DB into memory..."
    s = """ SELECT diaSourceId, opSimId, ssmId, 
            ra, decl, taiMidPoint, mag, snr  FROM
            %s.%s; """ % (mopsDatabases.DIAS_DB, mopsDatabases.DIAS_TABLE)
    cursor.execute(s)
    results = cursor.fetchall()
    for row in results:
        toRet[row[0]] = row

    print "... Done fetching dias."
    return toRet


def writeHeader(outf):
    outf.write("#!trackNum diaSourceId opSimId ssmId ra decl taiMidPoint mag snr\n")



def writeTrackDias(trackId, diasData, outf):
    for dia in diasData:
        trackString = ""
        trackString += str(trackId) + " "
        for datum in dia:
            trackString += str(datum) + " "
        outf.write(trackString + "\n")
        

def tracksToLynneFormat(inf, outf, cursor):
    idToDataMap = fetchAllDiasFromDb(cursor)
    line = inf.readline()
    count = 0
    writeHeader(outf)
    while line != "":
        diaIds = map(int, line.split())
        diasData = getDiaInfo(diaIds, idToDataMap)
        writeTrackDias(count, diasData, outf)
        count += 1
        if count % 1000 == 0:
            print "Written ", count, " tracks so far."
        line = inf.readline()



if __name__=="__main__":
    inf = file(sys.argv[1],'r')
    outf = file(sys.argv[2],'w')
    curs = mopsDatabases.getCursor()

    tracksToLynneFormat(inf, outf, curs)

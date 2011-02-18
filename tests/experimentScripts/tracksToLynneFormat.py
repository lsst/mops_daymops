#!/usr/bin/env python

""" jmyers feb 18 2011

convert tracks (or tracklets even) from byDiaId format (one line per
track, white-space delimited DiaIds which build the track) into
something like the format of Lynne's original DiaSource dump:

#!trackNum diaSourceId opSimId ssmId ra decl taiMidPoint mag snr

Tracks will be assigned IDs sequentially starting at zero.
"""



import mopsDatabases
import sys
import MySQLdb as db


def getDiaInfo(diaIds, cursor):
    s = """ SELECT diaSourceId,  opSimId, ssmId, ra, decl, taiMidPoint,
    mag, snr FROM %s.%s WHERE diaSourceID IN ( """ % (
        mopsDatabases.DIAS_DB, mopsDatabases.DIAS_TABLE)
    first = True
    for d in diaIds:
        if not first:
            s += ", "
        s += "%d" % d
        first = False
    s += ");" 
    cursor.execute(s)
    return cursor.fetchall()


def writeHeader(outf):
    outf.write("#!trackNum diaSourceId opSimId ssmId ra decl taiMidPoint mag snr")



def writeTrackDias(trackId, diasData, outf):
    trackString = ""
    for dia in diasData:
        trackString += str(trackId) + " "
        for datum in dia:
            trackString += str(datum) + " "
        outf.write(trackString)

def tracksToLynneFormat(inf, outf, cursor):
    line = inf.readline()
    count = 0
    writeHeader(outf)
    while line != "":
        diaIds = map(int, line.split())
        diasData = getDiaInfo(diaIds, cursor)
        writeTrackDias(count, diasData, outf)
        count += 1
        line = inf.readline()



if __name__=="__main__":
    inf = file(sys.argv[1],'r')
    outf = file(sys.argv[2],'w')
    curs = mopsDatabases.getCursor()

    tracksToLynneFormat(inf, outf, curs)

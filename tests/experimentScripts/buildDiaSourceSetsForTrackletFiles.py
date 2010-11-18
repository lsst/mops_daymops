#!/usr/bin/env python

"""

jon myers oct 27

given tracklets (or track) files expressed as sets of DiaSourceIds,
build files containing ONLY THE DIAS REFERENCED IN THE GIVEN FILE.

Currently uses flat-file format like Lynne's dump of DiaSources in the format 

"""


import MySQLdb as db

DB_USER="jmyers"
DB_PASSWD="jmyers"
DB_HOST="localhost"

DIAS_DB="mops_noDeepAstromError"
DIAS_TABLE="fullerDiaSource"



def lookUpDia(diaId, cursor=None):
    s = """ SELECT diaSourceId, opSimId, ssmId, ra, decl, taiMidPoint, mag, snr FROM 
                %s.%s WHERE diaSourceId=%d; """ % (DIAS_DB, DIAS_TABLE, diaId)
    cursor.execute(s)
    [diaId, obsHistId, ssmId, ra, dec, mjd, mag, snr] = cursor.fetchone()
    d = DiaSource(diaId=diaId, obsHistId=obsHistId, ssmId=ssmId, ra=ra, dec=dec, mjd=mjd, mag=mag, snr=snr)
    return d





class DiaSource:
    def __init__(self, diaId, obsHistId, mjd, ssmId, ra, dec, mag, snr):
        self.diaId = diaId
        self.obsHistId = obsHistId
        self.obsTime = mjd
        self.ssmId = ssmId
        self.ra = ra
        self.dec = dec
        self.mag = mag
        self.snr = snr

    def getDiaId(self):
        return self.diaId

    def getObsTime(self):
        return self.obsTime

    def getSsmId(self):
        return self.ssmId

    def getObsHistId(self):
        return self.obsHistId
    
    def getRa(self):
        return self.ra

    def getDec(self):
        return self.dec

    def getMag(self):
        return self.mag

    def getSnr(self):
        return self.snr




def getAllDiasReferencedInFile(f):
    line = f.readline()
    toRet = set()
    while line != "":
        dias = map(int, line.split())
        for dia in dias:
            toRet.add(dia)
        line = f.readline()
    return toRet




def diaToString(dia):
    """ format a dia to Lynne's database dump format:

    diaId obsHistId ssmId ra dec expMjd mag snr"""

    return "%d %d %r %3.10f %3.10f %3.10f %3.10f %3.10f\n" % \
        (dia.getDiaId(), dia.getObsHistId(), dia.getSsmId(), dia.getRa(), dia.getDec(), dia.getObsTime(), dia.getMag(), dia.getSnr())




def writeDiasToFile(diaIds, dbcurs, outfile):
    for diaId in diaIds:
        dia = lookUpDia(diaId, dbcurs)
        outfile.write(diaToString(dia))




if __name__=="__main__":

    import sys

    dbconn = db.connect(host=DB_HOST, user=DB_USER, passwd=DB_PASSWD)
    dbcurs = dbconn.cursor()

    if len(sys.argv) != 3:
        raise Exception("""USAGE: buildDiaSourceSetsForTrackletFiles <trackletsFile><outputDias>""")
    
    trackletsFile = file(sys.argv[1],'r')
    outputDiaSources = file(sys.argv[2],'w')

    print "Getting set of all diaSources referenced in file..."
    allDiasReferenced = getAllDiasReferencedInFile(trackletsFile)
    print "Found ", len(allDiasReferenced) , " dias to write to outfile."

    print "Writing outfile..."
    writeDiasToFile(allDiasReferenced, dbcurs, outputDiaSources)
    print "Done."

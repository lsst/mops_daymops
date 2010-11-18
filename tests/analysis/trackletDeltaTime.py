#!/usr/bin/env python


"""

jmyers nov 12 

takes a set of tracklets files (in the by-diaId format) and a
connection to the database and calculates the time separation
end-to-end of the tracklets.

"""



import sys


def getDt(ids, allDias):
    dets = [allDias[i] for i in ids]
    mjds = map(lambda x: x.getObsTime(), dets)
    return max(mjds) - min(mjds)

def getDtsToOutfile(inTrackletsFile, outDtsFile, allDias):
    line = inTrackletsFile.readline()
    while line != "":
        ids = map(int, line.split())
        dt = getDt(ids, allDias)
        outDtsFile.write("%2.8f\n" % (dt))
        line = inTrackletsFile.readline()




class DiaSource:
    def __init__(self, diaId=None, obsTime=None, ssmId=None, obsHistId=None, ra=None, dec=None, mag=None):
        self.diaId = diaId
        self.obsTime = obsTime
        self.ssmId = ssmId
        self.obsHistId = obsHistId
        self.ra = ra
        self.dec = dec
        self.mag = mag

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


def readDias(diasDataFile):

    """ reads a dump of dias, which include diaId expMjd ssmId
    obsHistId for every diaSource.  Returns a dict mapping diaId to
    data."""
    idToDias = {}
    line = diasDataFile.readline()
    while line != "":
        [diaId, mjd, ra, dec, mag, obscode, ssmId, elong, angle] = line.split()
        diaId = int(diaId)
        [ra, dec, mjd, mag, elong, angle] = map(float, [ra, dec, mjd, mag, elong, angle])
        
        idToDias[diaId] = DiaSource(diaId=diaId, obsTime=mjd, ssmId=ssmId, 
                                    ra=ra, dec=dec, mag=mag)
        
        line = diasDataFile.readline()

    return idToDias



if __name__=="__main__":
    

    inTracklets = file(sys.argv[1],'r')
    inMitifile = file(sys.argv[2], 'r')
    print "Reading dias from file ", inMitifile
    allDias = readDias(inMitifile)
    print "Done."
    outDts = file(sys.argv[3], 'w')
    print "Writing tracks deltaTimes to ", outDts
    getDtsToOutfile(inTracklets, outDts, allDias)
    print "Done."

#!/usr/bin/env python

"""

jmyers may 11 2010

Split up the "simplified" DIAsources by night and put them in separate
files.  These are assumed to have a temporary format of

ID obshistid SSM_ID RA DECL MJD MAG sn

as created by running < dias_nodeep_pt1 awk '{print $2, $6, $10, $13,
$36, $39, $59}' > dias.short and then using addUniqueIds.py to add an
initial ID column.

lower-case items are ignored.  

"""


NIGHT_START_OFFSET=-.7

import sys

def getNightNum(MJD):

    for guess in range(int(MJD) - 2, int(MJD) + 3):
        if MJD > guess + NIGHT_START_OFFSET and MJD < guess + 1 + NIGHT_START_OFFSET:
            return guess 
    raise Exception("Couldn't get a good guess of night num for MJD %d" % (MJD))
    


if __name__=="__main__":
    inf = file(sys.argv[1], 'r')

    line = inf.readline();
    while line != "":
        diaId, obshistId, ssmId, ra, decl, MJD, mag, snr = line.split()
        diaId, ssmId = map(int, [diaId, ssmId])
        ra, decl, MJD, mag = map(float, [ra, decl, MJD, mag])

        mitiLine = "%i %f %f %f %f 566 %i 0. 0." % (diaId, MJD, ra, decl, mag, ssmId)

        nightNum = getNightNum(MJD)
        outf = file(str(nightNum) + ".miti", 'aw')
        print>>outf, mitiLine

        line = inf.readline()

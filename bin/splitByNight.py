#!/usr/bin/env python

"""

jmyers may 11 2010

march 29 2012: Get rid of MITI format and write things in fullerDiaSource format.

Split up the "fullerDiaSource"-format DIAsources by night and put them in separate
files.  

apr. 3 2012: Also, write out a per-obsHist file which holds all dias from a given image.

"""

# note that midnight at LSST is MJD+0.125 (or 0.166) days
#  (MJD = integer at midnight UTC, Chile/LSST local time is 
#    UTC -4 hours in standard,-3 hours daylight savings. 
#    which translates to midnight @ LSST = MJD + 0.125/0.166)
# which means ~NOON to NOON observations cover 
#   ~'night'-0.35 to ~'night' + 0.65   (where 'night' = int(MJD) at midnight)
#    a gap would be okay because we're don't observe that close to noon

night_start = -0.35
night_end = 0.65

OBSCODE='807'

import sys
import os.path

def getNightNum(mjd):
    """Determine night number for any MJD."""
    night = night = int(mjd+0.5-0.12)    
    return night
    

if __name__=="__main__":


    if len(sys.argv)<2:
        print "Usage: splitByNight.py filename nightlyOutputDir byObsHistOutputDir"
        print "  where filename = the input diasource file "
        print "  dia sources broken up by night will go in nightlyOutputDir"
        print "  dia sources broken up by image will go in byObsHistOutputDir"
        sys.exit(1)

    infile = open(sys.argv[1], 'r')
    outDir1 = sys.argv[2]
    outDir2 = sys.argv[3]

    prev_night = None

    # Read diasources from input file.
    for line in infile:
        diaId, obshistId, ssmId, ra, decl, MJD, mag, snr = line.split()
        diaId, obshistId, ssmId = map(int, [diaId, obshistId, ssmId])
        ra, decl, MJD, mag = map(float, [ra, decl, MJD, mag])

        # Determine the night number of this particular diasource and write to that file.
        nightNum = getNightNum(MJD)
        # Open new output file if needed.
        if nightNum != prev_night:
            outfile = open(os.path.join(outDir1, str(nightNum) + ".dias"), "aw")
            prev_night = nightNum
        # Write output line.
        print>>outfile, line.rstrip()

        # now write to a by-obshist dir
        outfile2 = open(os.path.join(outDir2, str(obshistId) + ".dias"), "aw")
        print>>outfile2, line.rstrip()

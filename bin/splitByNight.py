#!/usr/bin/env python

"""

jmyers may 11 2010

march 29 2012: no longer write output in MITI, get rid of MITI forever!

Split up the "fullerDiaSource"-format DIAsources by night and put them in separate
files.  

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

def getNightNum(mjd):
    """Determine night number for any MJD."""
    night = night = int(mjd+0.5-0.12)    
    return night
    

if __name__=="__main__":


    if len(sys.argv)<2:
        print "Usage: splitByNight.py filename <[start night] [end night]>"
        print "  where filename = the input diasource file "
        print "  and start/end night are optional"
        print "   (but equal to the MJD nightNum of the interval you want to extract)."
        sys.exit()

    infile = open(sys.argv[1], 'r')
    
    start_interval = -99
    stop_interval = 1e10
    if len(sys.argv) > 2:
        start_interval = int(sys.argv[2])
        stop_interval = int(sys.argv[3])

    prev_night = -99

    # Read diasources from input file.
    for line in infile:
        diaId, obshistId, ssmId, ra, decl, MJD, mag, snr = line.split()
        diaId, ssmId = map(int, [diaId, ssmId])
        ra, decl, MJD, mag = map(float, [ra, decl, MJD, mag])
        # Determine the night number of this particular diasource.
        nightNum = getNightNum(MJD)
        if (nightNum >= start_interval) & (nightNum <= stop_interval):
            # Open new output file if needed.
            if nightNum != prev_night:
                outfile = open(str(nightNum) + ".dias", "aw")
                outfile2 = open(str(nightNum) + ".miti", "aw")
                prev_night = nightNum
            # Write output line.
            print>>outfile, line.rstrip()
            print>>outfile2, "%d %f %f %f %f %s %d 0.0 0.0 " % (diaId, MJD, ra, decl, mag, OBSCODE, ssmId)

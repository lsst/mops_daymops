#!/usr/bin/env python

""" make runScripts for linkTracklets. """

import glob
import os

INFILES_GLOB=glob.glob("../*.miti")
START_T_RANGE_FOR_INFILE=lambda inf: inf[:-4] + "start_t_range"
BASENAME_FOR_INFILE=lambda inf: os.path.basename(inf[:-5])
OUTF_FOR_INFILE=lambda inf: os.path.basename(inf[:-4] + "cmd.sh")


def writeRunScript(infile, startTRangeFile):
    outf = file(OUTF_FOR_INFILE(infile),'w')
    outS = """#!/usr/bin/bash

# a quick test to determine the approximate runtime/num false tracklets from a run with min_sup == 2 and 5 day windows

BN=""" + BASENAME_FOR_INFILE(infile) + """

CMD="$AUTON_DIR/linkTracklets_modified/linkTracklets_modified file ../$BN.miti indicesfile $BN.c.tracks.byIndices start_t_range `cat ../$BN.start_t_range`   acc_r 0.02 acc_d 0.02 fit_thresh  0.000000250000 min_sup 3 min_obs 6 plate_width .00000001"

echo Running: $CMD

/usr/bin/time -o $BN.c.linkTracklets.runtime $CMD | tee $BN.c.linkTracklets.runlog
"""
    outf.write(outS)
    outf.close()

if __name__=="__main__":
    for infile in INFILES_GLOB:
        print infile
        startTRangeFile = START_T_RANGE_FOR_INFILE(infile)
        print startTRangeFile

        writeRunScript(infile, startTRangeFile)

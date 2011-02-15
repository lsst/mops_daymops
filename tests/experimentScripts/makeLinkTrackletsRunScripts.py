#!/usr/bin/env python

""" make runScripts for linkTracklets. """

import glob
import os

INFILES_GLOB=glob.glob("../*.miti")
START_T_RANGE_FOR_INFILE=lambda inf: inf[:-4] + "start_t_range"
CPP_START_T_RANGE_FOR_INFILE=lambda inf: inf[:-4] + "date.start_t_range"
CPP_DETS_FOR_INFILE=lambda inf: inf[:-4] + ".dets"
CPP_IDS_FOR_INFILE=lambda inf: inf[:-4] + ".ids"
CPP_TRACKS_FOR_INFILE=lambda inf: inf[:-4] + ".cpp.tracks"
CPP_OUTF_FOR_INFILE=lambda inf: os.path.basename(inf[:-4] + ".cpp.cmd.sh")

BASENAME_FOR_INFILE=lambda inf: os.path.basename(inf[:-5])
OUTF_FOR_INFILE=lambda inf: os.path.basename(inf[:-4] + "c.cmd.sh")


def writeRunScript(infile, startTRangeFile):
    outf = file(OUTF_FOR_INFILE(infile),'w')
    outS = """#!/usr/bin/bash


BN=""" + BASENAME_FOR_INFILE(infile) + """

CMD="$AUTON_DIR/linkTracklets_modified/linkTracklets_modified file ../$BN.miti indicesfile $BN.c.tracks.byIndices start_t_range `cat ../$BN.start_t_range`   acc_r 0.02 acc_d 0.02 fit_thresh  0.000000250000 min_sup 3 min_obs 6 plate_width .00000001"

echo Running: $CMD

/usr/bin/time -o $BN.c.linkTracklets.runtime $CMD | tee $BN.c.linkTracklets.runlog
"""
    outf.write(outS)
    outf.close()


def writeCppRunScript(infile, startTRangeFile):
    outf = file(CPP_OUTF_FOR_INFILE(infile),'w')
    dets = CPP_DETS_FOR_INFILE(infile)
    ids = CPP_IDS_FOR_INFILE(infile)
    tracks = CPP_TRACKS_FOR_INFILE(infile)
    args = "-d " + dets + " -t " + ids + " -o " + tracks

    outS = """#!/usr/bin/bash


CMD="$MOPS_DAYMOPS_DIR/bin/linkTracklets """ + args + """

echo Running: $CMD

/usr/bin/time -o $BN.cpp.linkTracklets.runtime $CMD | tee $BN.cpp.linkTracklets.runlog
"""
    outf.write(outS)
    outf.close()


if __name__=="__main__":
    for infile in INFILES_GLOB:
        print infile
        startTRangeFile = START_T_RANGE_FOR_INFILE(infile)
        print startTRangeFile

        writeRunScript(infile, startTRangeFile)
        writeCppRunScript(infile, startTRangeFile)

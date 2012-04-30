#!/usr/bin/env python

""" make runScripts for linkTracklets. """

import glob
import os

INFILES_GLOB=glob.glob("../*.dets")
START_T_RANGE_FOR_INFILE=lambda inf: inf[:-4] + "start_t_range"
CPP_START_T_RANGE_FOR_INFILE=lambda inf: inf[:-4] + "date.start_t_range"
CPP_DETS_FOR_INFILE=lambda inf: inf[:-4] + "dets"
CPP_IDS_FOR_INFILE=lambda inf: inf[:-4] + "ids"
CPP_TRACKS_FOR_INFILE=lambda inf: os.path.basename(inf)[:-4] + "cpp.tracks"
CPP_OUTF_FOR_INFILE=lambda inf: os.path.basename(inf[:-4] + "cpp.cmd.sh")

BASENAME_FOR_INFILE=lambda inf: os.path.basename(inf[:-5])
OUTF_FOR_INFILE=lambda inf: os.path.basename(inf[:-4] + "c.cmd.sh")

VTREE_THRESH=.0004



def writeCppRunScript(infile, startTRangeFile, vtree_thresh):
    outf = file(CPP_OUTF_FOR_INFILE(infile),'w')
    dets = CPP_DETS_FOR_INFILE(infile)
    ids = CPP_IDS_FOR_INFILE(infile)
    tracks = CPP_TRACKS_FOR_INFILE(infile)
    startRange = CPP_START_T_RANGE_FOR_INFILE(infile)
    logFile = BASENAME_FOR_INFILE(infile) + ".cpp.runlog"

    args = "-d " + dets + " -t " + ids + " -o " + tracks + " -F `cat " + startRange + "`" + " -e " + str(vtree_thresh)

    outS = """#!/usr/bin/bash


CMD="$MOPS_DAYMOPS_DIR/bin/linkTracklets """ + args + """ "

echo Running: $CMD

$CMD | tee """ + logFile + """
"""
    outf.write(outS)
    outf.close()


if __name__=="__main__":
    for infile in INFILES_GLOB:
        print infile
        startTRangeFile = START_T_RANGE_FOR_INFILE(infile)
        print startTRangeFile
        
        writeCppRunScript(infile, startTRangeFile, VTREE_THRESH)

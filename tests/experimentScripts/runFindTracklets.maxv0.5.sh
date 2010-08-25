#!/bin/bash


FINDTRACKLETS=$MOPS_DAYMOPS_DIR/bin/findTracklets

NIGHTLY_DIASOURCES="/workspace1/jmyers/nightlyDiasAstromErr/*.miti"
OUTPUT_DIR="/workspace1/jmyers/nightlyDiasAstromErr/tracklets/"

for NIGHTLY in $NIGHTLY_DIASOURCES
do
    CMD="$FINDTRACKLETS -i $NIGHTLY -o $OUTPUT_DIR/`basename $NIGHTLY .miti`.maxv0.5.tracklets -v .5 -m 0.0"
    echo running $CMD
    /usr/bin/time -o $OUTPUT_DIR/`basename $NIGHTLY .miti`.maxv0.5.time $CMD
    echo ""
    echo ""
done
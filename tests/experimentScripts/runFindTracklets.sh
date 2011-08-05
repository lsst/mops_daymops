#!/bin/bash


FINDTRACKLETS=$MOPS_DAYMOPS_DIR/bin/findTracklets

NIGHTLY_DIASOURCES=$PWD/*miti
OUTPUT_DIR=$PWD/tracklets 

MAXV=$1

if [ "$MAXV" == ""]
then
    print "USAGE: runFindTracklets.sh (max velocity)" 
    exit 1
fi

echo "Placing output data in $OUTPUT_DIR"

for NIGHTLY in $NIGHTLY_DIASOURCES
do
    CMD="$FINDTRACKLETS -i $NIGHTLY -o $OUTPUT_DIR/`basename $NIGHTLY .miti`.maxv0.5.tracklets -v ${MAXV} -m 0.0"
    echo running $CMD
    /usr/bin/time -o $OUTPUT_DIR/`basename $NIGHTLY .miti`.maxv0.5.time $CMD
    echo ""
    echo ""
done
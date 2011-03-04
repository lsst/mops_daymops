#!/bin/bash


FINDTRACKLETS=$MOPS_DAYMOPS_DIR/bin/findTracklets

NIGHTLY_DIASOURCES=../*.miti
OUTPUT_DIR=$PWD

if [ "$1" == "" ]
then
    echo "USAGE: runFindTracklets.sh <maxV>"
    exit 1
fi

MAXV=$1

echo "Placing output data in $OUTPUT_DIR"

for NIGHTLY in $NIGHTLY_DIASOURCES
do
    CMD="$FINDTRACKLETS -i $NIGHTLY -o $OUTPUT_DIR/`basename $NIGHTLY .miti`.maxv${MAXV}.tracklets -v ${MAXV} -m 0.0"
    echo running $CMD
    /usr/bin/time -o $OUTPUT_DIR/`basename $NIGHTLY .miti`.maxv{$MAXV}.time $CMD
    echo ""
    echo ""
done
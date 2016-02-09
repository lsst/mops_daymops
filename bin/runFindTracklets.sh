#!/bin/bash


FINDTRACKLETS=$MOPS_DAYMOPS_DIR/bin/findTracklets

echo "Working directory: $PWD"
NIGHTLY_DIASOURCES=`ls $PWD | grep .dias`
echo "Got diaSource detections: $NIGHTLY_DIASOURCES"
OUTPUT_DIR=$PWD/

MAXV=$1

if [ "$MAXV" == "" ]
then
    echo "USAGE: runFindTracklets.sh (max velocity)" 
    exit 1
fi

echo "Placing output data in $OUTPUT_DIR"

for NIGHTLY in $NIGHTLY_DIASOURCES
do
    CMD="$FINDTRACKLETS -i $NIGHTLY -o $OUTPUT_DIR/`basename $NIGHTLY .dias`.tracklets -v ${MAXV} -m 0.0"
    echo running $CMD
    LOGFILE=$OUTPUT_DIR/`basename $NIGHTLY .miti`.findTracklets.log
    $CMD  | tee $LOGFILE
    echo ""
    echo ""
done
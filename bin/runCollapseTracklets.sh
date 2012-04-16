#!/bin/bash

set -e

MOPS_BIN=$MOPS_DAYMOPS_DIR/bin/

for IN_TRACKLETS in *.tracklets
do
    echo "Processing $IN_TRACKLETS, starting at:"
    date
    BN=`basename $IN_TRACKLETS .tracklets`
    DIASFILE=../`basename $IN_TRACKLETS .tracklets`.dias
    python $MOPS_BIN/idsToIndices.py $IN_TRACKLETS $DIASFILE $IN_TRACKLETS.byIndices
    
    $MOPS_DAYMOPS_DIR/bin/collapseTracklets $DIASFILE $IN_TRACKLETS.byIndices \
	.002 .002 5 .05 $IN_TRACKLETS.collapsed | tee $BN.collapseTracklets.log

    $MOPS_DAYMOPS_DIR/bin/purifyTracklets \
	--detsFile $DIASFILE \
	--pairsFile $IN_TRACKLETS.collapsed \
	--outFile $IN_TRACKLETS.collapsed.pure | tee $BN.purifyTracklets.log

    $MOPS_DAYMOPS_DIR/bin/removeSubsets \
	--inFile $IN_TRACKLETS.collapsed.pure \
	--outFile $IN_TRACKLETS.final.byIndices  | tee $BN.removeSubsets.log

    python $MOPS_BIN/indicesToIds.py $IN_TRACKLETS.final.byIndices \
	$DIASFILE $IN_TRACKLETS.final.byDiaIds

done
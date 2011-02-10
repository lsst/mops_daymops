#!/bin/bash

set -e

MOPS_HACKS=$MOPS_DAYMOPS_DIR/tests/experimentScripts/

for IN_TRACKLETS in *.tracklets
do
    echo "Processing $IN_TRACKLETS, starting at:"
    date

    MITIFILE=../`basename $IN_TRACKLETS .maxv0.5.tracklets`.miti
    python $MOPS_HACKS/idsToIndices.py $IN_TRACKLETS $MITIFILE $IN_TRACKLETS.byIndices

    /usr/bin/time -o $BN.collapse.time \
	$MOPS_DAYMOPS_DIR/bin/collapseTracklets $MITIFILE $IN_TRACKLETS.byIndices \
	.002 .002 5 .05 $IN_TRACKLETS.collapsed

    /usr/bin/time -o $BN.purify.time \
	$MOPS_DAYMOPS_DIR/bin/purifyTracklets \
	--detsFile $MITIFILE \
	--pairsFile $IN_TRACKLETS.collapsed \
	--outFile $IN_TRACKLETS.collapsed.pure

    /usr/bin/time -o $BN.removeSubs.time \
	$MOPS_DAYMOPS_DIR/bin/removeSubsets \
	--inFile $IN_TRACKLETS.collapsed.pure \
	--outFile $IN_TRACKLETS.final.byIndices

    python $MOPS_HACKS/indicesToIds.py $IN_TRACKLETS.final.byIndices \
	$MITIFILE $IN_TRACKLETS.final.byDiaIds

done
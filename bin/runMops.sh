#!/bin/bash

# jmyers - a quick script that should take a set of diaSources, a name
# for the DB an DB table, and go ahead and populate the DB with
# diaSources and run find/collapse/linkTracklets on the data.

set -e

DIAS_FILE=`readlink -f "${1}" `
#DIAS_FILE="$1"


#set these if you care to
WINDOW_SIZE=15
TRACKLET_MAXV=.5
DB_HOST="localhost"
DB_USER="jmyers"
DB_PASS="jmyers"


echo "Working on dias from " $DIAS_FILE  | tee -a mopsSetup.log
echo "Using MOPS software from " $MOPS_DAYMOPS_DIR | tee -a mopsSetup.log

if [ "$DIAS_FILE" = "" ] || [ "$MOPS_DAYMOPS_DIR" = "" ]
then
    echo "USAGE: runMops.sh <input dias file, Lynne format> "
    echo "you MUST have mops_daymops set up by eups first."
    exit 1
fi


echo " EUPS PACKAGES: " | tee -a mopsSetup.log
eups list | tee -a mopsSetup.log

# now do the real work.

# split up the giant input file into manageable chunks.
mkdir diasByObsHist/
python $MOPS_DAYMOPS_DIR/bin/splitByNight.py $DIAS_FILE ./ diasByObsHist/


#now find the tracklets.
mkdir tracklets
cd tracklets

bash $MOPS_DAYMOPS_DIR/bin/runFindTracklets.sh $TRACKLET_MAXV
# collapseTracklets
bash $MOPS_DAYMOPS_DIR/bin/runCollapseTracklets.sh

# bin tracklets in prep for building linkTracklets infiles
mkdir byObsHistId/

for TRACKLETS in *.tracklets.final.byDiaIds
do
    DIASFILE=../`basename $TRACKLETS .tracklets.final.byDiaIds`.dias
    python $MOPS_DAYMOPS_DIR/bin/binTrackletsByStartImage_noDb.py $TRACKLETS $DIASFILE byObsHistId/ 
done

exit 0

mkdir linkTrackletsInfiles

python $MOPS_DAYMOPS_DIR/bin/makeLinkTrackletsInput_byImages.py \
 $PWD/byObsHistId/ $PWD/linkTrackletsInfiles/ $DIAS_DB $DIAS_TABLE $WINDOW_SIZE

cd linkTrackletsInfiles
mkdir tracks
cd tracks
python $MOPS_DAYMOPS_DIR/bin/makeLinkTrackletsRunScripts.py

for CMD in *.cpp.cmd.sh
do
   bash $CMD
done

# grade the tracks we got out 
python $MOPS_DAYMOPS_DIR/bin/../analysis/countTrueFalseTracksByObsIdAndFindableObjects.py $DIAS_FILE \*.tracks 


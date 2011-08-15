#!/bin/bash

# jmyers - a quick script that should take a set of diaSources, a name
# for the DB an DB table, and go ahead and populate the DB with
# diaSources and run find/collapse/linkTracklets on the data.

set -e

DIAS_FILE=`readlink -f "${1}" `
#DIAS_FILE="$1"
# we expect DIAS_DB NOT exist already.
DIAS_DB="$2"
# DIAS_TABLE should NOT exist already.
DIAS_TABLE="$3"


#set these if you care to
WINDOW_SIZE=30
TRACKLET_MAXV=.75
DB_HOST="localhost"
DB_USER="jmyers"
DB_PASS="jmyers"


echo $DIAS_FILE $DIAS_DB $DIAS_TABLE $MOPS_DAYMOPS_DIR
if [ "$DIAS_FILE" == ""] || [ "$DIAS_DB" == "" ] || [ "$DIAS_TABLE" == "" ] || [ "$MOPS_DAYMOPS_DIR" == "" ]
then
    echo "USAGE: runMops.sh <input dias file, Lynne format> <dias DB, should NOT exist already> <dias Table to create (should NOT exist already.)>"
    echo "you MUST have mops_daymops set up by eups first."
    exit 1
fi


#write some helpful data

# for some reason, svn info $MOPS_DAYMOPS_DIR doesn't work but cd-ing and running svn info does.
OLDPWD=$PWD
cd $MOPS_DAYMOPS_DIR
echo " CURRENT SVN VERSION: " | tee -a $OLDPWD/mopsSetup.log
svn info  | tee -a $OLDPWD/mopsSetup.log
cd $OLDPWD

echo "" | tee -a mopsSetup.log
#echo "Creating DB: $DIAS_DB " | tee -a mopsSetup.log
#echo "Reading DiaSources from : $DIAS_FILE " | tee -a mopsSetup.log
#echo "Copying DiaSources to table : $DIAS_TABLE " | tee -a mopsSetup.log
echo "" | tee -a mopsSetup.log
echo " EUPS PACKAGES: " | tee -a mopsSetup.log
eups list | tee -a mopsSetup.log

# now do the real work.

MOPS_HACKS="$MOPS_DAYMOPS_DIR/tests/experimentScripts/"


# set up the databases for later.
SQL="mysql -h $DB_HOST -u $DB_USER -p$DB_PASS"

#echo "CREATE DATABASE $DIAS_DB; USE $DIAS_DB; `cat $MOPS_HACKS/fullerDiaSource.sql`;" | $SQL
#echo "CREATE TABLE $DIAS_TABLE LIKE fullerDiaSource;" | $SQL $DIAS_DB
#echo "LOAD DATA INFILE '$DIAS_FILE' INTO TABLE 
#      $DIAS_TABLE FIELDS TERMINATED BY ' '" | $SQL $DIAS_DB

python  $MOPS_HACKS/splitByNight.py $DIAS_FILE

mkdir tracklets

bash $MOPS_HACKS/runFindTracklets.sh $TRACKLET_MAXV

cd tracklets

# collapseTracklets
bash $MOPS_HACKS/runCollapseTracklets.sh

# bin tracklets in prep for building linkTracklets infiles
mkdir byObsHistId/

for TRACKLETS in *.tracklets.final.byDiaIds
do
    python $MOPS_HACKS/binTrackletsByStartImage.py $TRACKLETS byObsHistId/ $DIAS_DB $DIAS_TABLE
done

mkdir linkTrackletsInfiles

python $MOPS_HACKS/makeLinkTrackletsInput_byImages.py \
 $PWD/byObsHistId/ $PWD/linkTrackletsInfiles/ $DIAS_DB $DIAS_TABLE $WINDOW_SIZE

cd linkTrackletsInfiles
mkdir tracks
cd tracks
python $MOPS_HACKS/makeLinkTrackletsRunScripts.py

for CMD in *.cpp.cmd.sh
do
   bash $CMD
done

# grade the tracks we got out 
python $MOPS_HACKS/../analysis/countTrueFalseTracksByObsIdAndFindableObjects.py $DIAS_FILE \*.tracks 


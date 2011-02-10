Here's a greatly simplified guide to running MOPS at the moment, which
will run findTracklets, collapseTracklets linkTracklets for you.

I was using Bash when I came up with these, you may need to change a
few things if you're using *csh.

# set up your environment
setlsst
setup mysqlpython
setup mops_daymops
export MOPS_HACKS=$MOPS_DAYMOPS_DIR/tests/experimentScripts/

# get data
mkdir myMopsRun
cd myMopsRun
wget --user=USER --password=PASSWORD dias_pt1_nodeep.short.astromErr


# populate the DB for later. I assume you have the OpSim DB already.
echo "CREATE DATABASE myMops; USE myMops; `cat fullerDiaSource.sql`;" | mysql
echo "LOAD DATA INFILE '$PWD/dias_pt1_nodeep.short.astromErr' INTO TABLE 
      fullerDiaSource FIELDS TERMINATED BY ' '" | mysql myMops


# findTracklets
python  /home/jmyers/sandbox/daymops_trunk/tests/experimentScripts/splitByNight.py dias_pt1_nodeep.short.astromErr

mkdir tracklets

bash $MOPS_HACKS/runFindTracklets.maxv0.5.sh

cd tracklets

# collapseTracklets
bash $MOPS_HACKS/runCollapseTracklets.sh


##################################
# make linkTracklets input files #
##################################

# edit mopsDatabases.py at the top to make sure
# it points to the right DBs and tables for opSim and
# your new myMops DB, etc.

emacs $MOPS_HACKS/mopsDatabases.py

# bin the tracklets by their source images.

mkdir byObsHistId/

for TRACKLETS in *.tracklets.final.byDiaIds
do
    python $MOPS_HACKS/binTrackletsByStartImage.py $TRACKLETS byObsHistId/
done

mkdir linkTrackletsInfiles

# edit makeLinkTrackletsInput_byImages.py so that it knows where to
# look for your byObsHistId directory, and your linkTrackletsInfiles
# directory.  These are set by constants at the top of the file.
emacs $MOPS_HACKS/makeLinkTrackletsInfiles_byObsHist.py

python $MOPS_HACKS/makeLinkTrackletsInfiles_byObsHist.py

cd linkTrackletsInfiles


# now make the run scripts which execute linkTracklets

mkdir tracks
cd tracks
python makeLinkTrackletsRunScripts.py

#########################
# RUNNING LINKTRACKLETS #
#########################

for LINKTRACKLETS_RUN_CMD in *.cmd.sh
do
   bash $LINKTRACKLETS_RUN_CMD
done

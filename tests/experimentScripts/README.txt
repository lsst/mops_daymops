Jmyers Aug 23

The following is a set of instructions for running
find/collapse/linkTracklets on some diaSources, the hackish non-LSST
way.

All the scripts should be in the same directory as this readme file.




BUILDING C++ FIND/LINKTRACKLETS (etc.)
---------------------------------------

You'll need to first install the LSST build system and at least
pex_exceptions.  Instructions for installing the LSST build system and
the full LSST pipeline development tools are available at
http://dev.lsstcorp.org/trac/wiki/Installing

Then you'll want to declare DayMOPS to EUPS and build it using SCons:

# get the code from SVN if you haven't already
$ svn co svn+ssh://lsstarchive.ncsa.uiuc.edu/DMS/mops/daymops/trunk/ ~/sandbox/daymops_trunk

# declare mops_daymops (version trunk, -c for current)
$ eups declare -c -r ~/sandbox/daymops_trunk   mops_daymops trunk

# set up mops_daymops - this will make sure all dependencies are installed 
$ setup mops_daymops

# now go build with scons
$ cd ~/sandbox/daymops_trunk
$ scons   # or scons opt=3 for a faster version.




GETTING SOME DIASOURCE DATA
---------------------

Lynne has data available at UW:
http://www.astro.washington.edu/users/ljones/ephems/

I am using
http://www.astro.washington.edu/users/ljones/ephems/dias_pt1_nodeep.short.astromErr

which has the format: 

#opsimID   ssmid  ra  dec   expmjd  mag  SNR

We want the middle five. We will also need to assign unique IDs to
each DiaSource so use the addUniqueIds script:

$ python $MOPS_HACKS/addUniqueIds.py dias_pt1_nodeep.short.astromErr
    dias_pt1_nodeep.short.astromErr.plusIds







PREPARING FOR FINDTRACKLETS
-----------------------------

Next you'll need to split the data up by night.  Use
splitByNight.py. If your input file has a different format you'll need
to do some hacking.  The nightly files will be created in your current directory, so make sure that your working directory is set intelligently.

$ cd nightlyDiasAstromErr

$ python $MOPS_HACKS/splitByNight.py dias_pt1_nodeep.short.astromErr.plusIds

This will take about 5 or 10 minutes, maybe a bit more.

Output will be in PanSTARRS MITI format: 
ID EPOCH_MJD RA_DEG DEC_DEG MAG OBSCODE OBJECT_NAME LENGTH
ANGLE 

OBSCODE LENGTH and ANGLE are all set to dummy values as they are not
used by our sky-plane algorithms.  They *will* be relevant when we get
to IOD though so beware!





RUNNING FINDTRACKLETS
----------------------------

Make sure you have dayMops set up via Eups.
$ source /path/to/your/loadLSST.sh
$ setup mops_daymops

Modify the first three variables in runFindTracklets.maxv0.5.sh to
suit your needs.  This script assumes you're using the C++
findTracklets.  Currently, C++ findTracklets takes a limited set of
command-line arguments; to adjust the fine-grained functionality
you'll have to modify findTracklets.h and recompile (there should be
documentation inside of findTracklets.h)

The script will also measure timing information for you.

$ bash $MOPS_HACKS/runFindTracklets.maxv0.5.sh

This will take a bit of time - perhaps an hour.

Currently C++ findTracklets writes its output as ASCII text, one
tracklet per line, with tracklets expressed as the set of Dia IDs of the diaSources in the tracklet, e.g.

----start example---------
1 2
4 10
5 7
-----end example---------


C findTracklets uses a different output format.





OPTIONAL PHASE: RUNNING COLLAPSETRACKLETS
--------------------------------------

TBD.




PREPARING FOR LINKTRACKLETS
-----------------------------

LinkTracklets uses bundled up sets of tracklets covering ~15 days of
observations and looks for tracks which start and end sometime within
that 15 days.  Depending on which flavor (C Auton / C++ LSST) of
linkTracklets you use the input formats are different.

Use make15DayWindows.py to bundle up your tracklets for C++
linkTracklets.  There is another script which will convert the C++
input files to work with C linkTracklets if you want to go that route.

Configuration of make15DayWindows is done by setting the constants at
the top of the file.  Change the constants so they point to your
diaSource files (.miti) and your tracklets.   Make sure that you create the output directory if needed:

e.g.
$ mkdir /workspace1/jmyers/nightlyDiasAstromErr/15DayWindowsMaxv0.5/

Then run the script:

$ python $MOPS_HACKS/make15DayWindows.py

Go get yourself a coffee, this will take ~20 min.

When you get your outfiles, they will be .dets files (sets of
diaSources in the window) and .ids files, which are ASCII text,
newline-delimited sets of LINE NUMBERS into the .dets file describing
each tracklet.*

The file names should be like:

night_49544_through_49557.dets
night_49544_through_49557.ids



*: LinkTracklets essentially wants to build a big array of detections
and an represents tracklets as an array of sets of indexes into the
detections array.  C linkTracklets does this invisibly, but to speed
things up C++ linkTracklets expects an external tool like
make15DayWindows.py to do this - hence we represent tracklets as sets
of line numbers.



RUNNING C++ LINKTRACKLETS
-----------------------------


$ $MOPS_DAYMOPS_DIR/bin/linkTracklets \
   -d night_49544_through_49557.dets \
   -t night_49544_through_49557.ids \
   -o night_49544_through_49557.tracks

For now, more advanced settings can be configured by editing linkTracklets.h and recompiling.

At the moment, C++ linkTracklets will look for ALL tracks which appear
within the entire window of time.  This can lead to some redundant
tracks - for example, a track which starts on night 49641 and ends on
night 49644 will be found while searching on the windows
night_49640_through_49654 as well as night_49641_through_49654.

Adding the option to limit searching to tracks starting or ending on a
subset of the nights in the data set is a forthcoming improvement
which should come Real Soon Now.







PREPARING FOR C LINKTRACKLETS
------------------------------

If you really must...  Then here are instructions.  It's going to get
awfully ugly, so take a deep breath and try not to get lost. 

To use the original Auton C implementation of linkTracklets, you'll
need to convert input file formats.  For the time being, C
linkTracklets supports searching for tracks which start or end on a
subset of the input data.

You'll want to make sure you're using the modified linkTracklets which
writes its output continuously, otherwise your linkTracklets will
probably crash.

C linkTracklets takes input as a single MITI file.  Use
trackletsToMiti to convert one .dets and .ids file (as would be input
to linkTracklets) to a single giant MITI file, e.g.:

$ python $MOPS_HACKS/trackletsToMiti.py night_49544_through_49557.dets \
   night_49544_through_49557.ids \
   night_49544_through_49557.tracklets.miti

You may wish to convert all of your files this way - you might want to
try using a BASH script like this:

for DETS in *.dets 
do
   BN=`basename $DETS .dets`
   python $MOPS_HACKS/trackletsToMiti.py $DETS $BN.ids $BN.tracklets.miti
done

Then go for a jog or find some good reading material.  The MITI format
generates some nasty bloat so there will be a lot of disk activity.
Which is funny, because C linkTracklets has to go back at runtime and
convert its in-memory representation of the data *back* to a format
more like the original one we were using as input to C++
linkTracklets anyway.

Note that the input file format for C linkTracklets is "lossy": it
will represent input as one giant MITI file, with the ID column used
for Tracklet IDs.  Contiguous detections in the MITI file which share
the same ID are trusted to be in the same tracklet.  The script
makeCheatSheetForMITI.py can be used to generate a "cheat sheet" file
such that line X of the "cheat sheet" file contains the DiaSourceId of
the DiaSource at line X of the .miti.tracklets file.  Confused yet?
See "interpreting the output of C linktracklets" for more info.






RUNNING C LINKTRACKLETS
--------------------------------

Anyway, once you've used trackletsToMiti.py to get your
.tracklets.miti files, you can run C linkTracklets:

$ AUTON_DIR/linkTracklets/linkTracklets \
  file night_49557_through_49572.tracklets.miti \
  indicesfile night_49557_through_49572.c.tracks.byIndices  \
  acc_r 0.02 acc_d 0.02 

acc_r and acc_d are the maximum accelerations of tracks for which you
will search, specified in RA and Dec.


There is one big upside to C linkTracklets right now: given a window,
it takes an optional argument to choose the set of track start and/or
end times for which it will attempt searching.  So if you have a data
set from nights 49563 through 49572 and another data set from nights
49564 through 49572, you can set linkTracklets to search *only* for
tracks starting on night 49563 when looking at the first data set.

Use the script makeRunScripts.py to generate BASH shell scripts which
will run linkTracklets while avoiding redundant searching.  It will do
this by peeking at the .tracklets.miti files and finding out when
their first tracklet time is.  Set the constants at the top of the
file to point to your .tracklets.miti files and then:

$ python $MOPS_HACKS/makeRunScripts.py





INTERPRETING THE OUTPUT OF C LINKTRACKLETS
----------------------------------------------

The "indicesfile" option is only supported by the modified version of
linkTracklets which I hacked together.  It will dump tracks
represented as ASCII text, tracks separated by newlines, and each
track will be represented as *A SET OF LINE NUMBERS* from the
.tracklets.miti file.  This is pretty ugly, but since I don't know my
way around Kubica's code that well (and because the input file format
lacks DiaSourceIds) it was the best thing we could come up with.

If you'd like to get the output in a more sane format, then you'll
need to first run the makeCheatSheetForMITI.py to get the DiaIds
corresponding to each entry in the MITI file, and you'll then need to
run indicesToIds.py to get tracks represented as sets of DiaSource
Ids.

Jmyers Aug 23

The following is a set of instructions for running
find/collapse/linkTracklets on some diaSources. 

In the future these scripts (or more likely, better versions of
all of this) will be modified so that pipelines can run each
stage of find/collapse/linkTracklets on particular sets of data.

All the scripts should be in the same directory as this readme file.



INSTALLING/BUILDING C++ FIND/LINKTRACKLETS (etc.)
---------------------------------------

Instructions online at http://dev.lsstcorp.org/trac/wiki/MOPS/Installing_MOPS



GETTING SOME INPUT DATA
---------------------

See instructions online at http://dev.lsstcorp.org/trac/wiki/MOPS/Data

Tests below run using example data :
http://www.astro.washington.edu/users/ljones/ephems/dias_pt1_nodeep.short.astromErr

which has the format: 

#diasourceID opsimID   ssmid  ra  dec   expmjd  mag  SNR


If you end up with diasource data which does not have the diasourceID, you
can add it using the script : 

$ python $MOPS_HACKS/addUniqueIds.py dias_pt1_nodeep.short.astromErr
    dias_pt1_nodeep.short.astromErr.plusIds

but this is not necessary using the most recent input data. 




PREPARING FOR FINDTRACKLETS
-----------------------------

Next you'll need to split the data up by night.  Use
splitByNight.py - this assumes that your input format
is as above: 
#diasourceID opsimID ssmID ra dec expmjd mag SNR

This script will create individual files of diasources in MITI
format, for each night. As there will be many different output
files, it's best to run this script in a working directory you
intend to use for this run of findTracklets/linkTracklets only.


$ cd nightlyDiasAstromErr

$ python $MOPS_HACKS/splitByNight.py dias_pt1_nodeep.short.astromErr

You can also choose to extract only certain nights from the entire range using

$ python $MOPS_HACKS/splitByNight.py dias_pt1_nodeep.short.astromErr 49639 49650 
(for example). 

Depending on the input data size, this will take a few minutes. 

Output will be in PanSTARRS MITI format: 
ID EPOCH_MJD RA_DEG DEC_DEG MAG OBSCODE OBJECT_NAME LENGTH
ANGLE 

LENGTH and ANGLE are all set to dummy values as they are not
used by our sky-plane algorithms. The OBSCODE is necessary for
doing orbit fitting (and has been set to approximately LSST positions). 



RUNNING FINDTRACKLETS
----------------------------

Make sure you have dayMops set up via Eups.
$ source /path/to/your/loadLSST.sh  (aka setlsst)
$ setup mops_daymops

You can now use $MOPS_HACKS/runFindTracklets.maxv.0.5.sh or 
$MOPS_HACKS/runFindTracklets.maxv.0.5.csh to run FindTracklets on all 
of the *miti files you just created. 

Modify the first three variables (defining the location of findTracklets 
and your input/output data) to suit your needs. For most uses of these 
scripts they'll be fine. 

This script assumes you're using the C++
findTracklets.  Currently, C++ findTracklets takes a limited set of
command-line arguments; to adjust the fine-grained functionality
you'll have to modify findTracklets.h and recompile (there should be
documentation inside of findTracklets.h)

The script will also measure timing information for you.

$ bash $MOPS_HACKS/runFindTracklets.maxv0.5.sh

This will take a bit of time - perhaps an hour. 

Note that you could restrict this to run ONLY on particular nights of
data by removing the other night's input from the source directory.

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

Go get yourself a coffee, this will take ~20 min (if you're working with
a lot of data). 

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





RUNNING C LINKTRACKLETS (the new way, distributed by images)
--------------------------------------

Our last set of runs on Abe demonstrated that a single night of
processing can require > 150 hours of CPU time. So to try to
distribute the workload differently, we now have the option of making
infiles with some maximum number of first (or perhaps someday last)
endpoint images in the files.

First you'll need to break up the tracklets by start image.

You WILL need the Opsim database as well as the a database of Dias
including the obsHistId of each image producing the diaSources.  

$ mkdir trackletsByStartImage
$ for TRACKLETSFILE in *.tracklets.byDiaId
do
   python $MOPS_HACKS/binTrackletsByStartImage.py trackletsByStartImage
done

You'll then need to make the new data sets for linkTracklets. They
will be automatically paired with the correct start_t_range in order
to limit the number of first images per data set.

Open makeLinkTrackletsInput_byImages.py and set all the constants at
the top to something sane for you; you'll mostly want to edit
TRACKLETS_BY_OBSHIST_DIR to be the location of trackletsByStartImage/
(as you just created above) and set up the database references
correctly.

Then update the MAX_START_IMAGES_PER_RUN settings and and
TRACKING_WINDOW_DAYS settings to your liking.


Then just run:

$ python $MOPS_HACKS/makeLinkTrackletsInput_byImages.py


This will create lots of input files with .start_t_range files.  To run C linkTracklets on data set foo:

$ $AUTON_DIR/linkTracklets_modified/linkTracklets file foo.miti start_t_range `cat foo.start_t_range` indicesfile foo.tracks.byIndices [other params - accel, RMS, etc.]




RUNNING C++ LINKTRACKLETS (distributed by night)
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







PREPARING FOR C LINKTRACKLETS (distributed by night)
------------------------------

If you really must...  Then here are instructions.  It's going to get
awfully ugly, so take a deep breath and try not to get lost. 
** However, we really would appreciate comparisons of the C and C++
linnkTracklets runs, so please do go ahead with this!

To use the original Auton C implementation of linkTracklets, you'll
need to convert input file formats.  For the time being, C
linkTracklets supports searching for tracks which start or end on a
subset of the input data.

You'll want to make sure you're using the modified linkTracklets which
writes its output continuously, from auton/linkTracklets_modified.

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






RUNNING C LINKTRACKLETS (distributed by night)
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

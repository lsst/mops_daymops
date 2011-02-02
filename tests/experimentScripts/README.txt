Jmyers Oct 22

Updated thoroughly to describe how I'm currently doing things.

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

You will also need the opsim database which was used to generate the
detections; currently this is opsim_3_61.




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

LinkTracklets uses bundled up sets of tracklets covering several days
of observations and looks for tracks which start and end sometime
within that window of time.  Depending on which flavor (C Auton / C++
LSST) of linkTracklets you use the input formats are different;
fortunately our new method generates both inputs simultaneously.

Use 


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





RUNNING LINKTRACKLETS (the new way, distributed by images)
--------------------------------------

Our last set of runs on Abe demonstrated that a single night of
processing can require > 150 hours of CPU time. So to try to
distribute the workload differently, we now have the option of making
infiles with some maximum number of first (or perhaps someday last)
endpoint images in the files.

You WILL need the Opsim database as well as the a database of Dias
including the obsHistId of each image producing the diaSources.  

First you'll need to break up the tracklets by start image.  This is
done with the script binTrackletsByStartImages.py; set the constants
at the top of the file to suit your needs. Then run it like  this:

Do that for every output file from findTracklets.  

$ mkdir trackletsByStartImage
$ for TRACKLETSFILE in *.tracklets.byDiaId
do
   python $MOPS_HACKS/binTrackletsByStartImage.py trackletsByStartImage
done


(NB: I think this only works with the output of C++ linkTracklets.)

There will be one file created for each ObsHistId of every image with
a tracklet "rooted" in that image.  All tracklets "rooted" in that
image will be grouped together.

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

$ $AUTON_DIR/linkTracklets_modified/linkTracklets file foo.miti start_t_range `cat foo.start_t_range` indicesfile foo.tracks.byIndices [args for other params - accel, RMS, etc.]







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



RUNNING IOD WITH ORBIT_SERVER (ORBFIT)
----------------------------------

Once you have your tracks in DiaId format, you can convert them to
input for orbit_server.x, Milani's IOD tool. 

Milani's IOD requires knowledge of per-detection astrometric error
bars. To calculate those, you'll need to use a database which holds
per-image information (e.g. the opsim table) and another table which
maps diaSources, including their parent images.

Use buildOrbitServerInput.py when you're ready; modify the
database/table names at the top of the file to suit your needs.  

orbit_server.x expects some maximum number of detections per image it
reads in input. This means that a large set of tracks may be broken up
into many input files for orbit_server.x.  Create a directory where
you'd like all those outfiles to go and cd there.

You'll need to tell the buildOrbitServerInput.py script a prefix for these outfiles, e.g. myTracks.  It will create many files, e.g.

myTracks_0.in.request
myTracks_0.in.tracklet
myTracks_0.in.manifest
myTracks_1.in.request
... etc.

So in the directory where you want all these files written: 

$ python $MOPS_HACKS/buildOrbitServerInput.py /path/to/tracks.byDiaId myTracks



Next, to run orbit_server on all these files, use runOrbitServer.sh.  Modify the file so it points to the orbit_server.x executable on your system, then just run it in the directory with the in.request,  in.tracklet, etc. files.


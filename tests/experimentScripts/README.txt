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

Next, you'll probably want to run collapseTracklets to turn your some
of your short, length-2 tracklets into longer tracklets.
CollapseTracklets really has several phases: first collapsing, then
purifying, then removing subset tracklets.

CollapseTracklets requires you to give some tolerances on the
command-line.  The ones I use are as follows: .002 .002 5 .05 (see
example commands below)

You should have a bunch of detections files like 49524.miti and
tracklets files with names like 49524.maxv0.5.tracklets, etc.  For
each one, you'll want to do something like this:

# first get the tracklets as a set of indices into the MITI file
$ python $MOPS_HACKS/idsToIndices.py 49524.maxv0.5.tracklets 49524.miti \
      49524.maxv0.5.tracklets.byIndices

# collapse the tracklets
$ $MOPS_DAYMOPS_DIR/bin/collapseTracklets 49524.miti 49524.maxv0.5.tracklets \
    .002 .002 5 .05 49524.maxv0.5.collapsed

# purify
$ $MOPS_DAYMOPS_DIR/bin/purifyTracklets --detsFile 49524.miti \
    --pairsFile 49524.maxv0.5.collapsed \
    --outFile 49524.maxv0.5.collapsed.pure

# remove subsets
$ $MOPS_DAYMOPS_DIR/bin/removeSubsets --inFile 49524.maxv0.5.collapsed.pure \
    --outFile 49524.maxv0.5.final.tracklets.byIndices 

# convert from by-file-index to by-dia-ids format
$ $MOPS_HACKS/indicesToIds 49524.maxv0.5.final.tracklets.byIndices \
  49524.miti 49524.maxv0.5.final.tracklets.byDiaIds

Since you have lots of nights, you'll probably want to use a shell
script or Python script to do the above on all your data (I used a
Bash for loop and the the 'basename' command.)



====================================
PREPARING FOR LINKTRACKLETS
====================================


REQUIRED ITEMS IN DATABASE:

In order to make the linkTracklets input and orbit_server input, some
scripts will expect the following things in the database:

1) the fullerDiaSource table, which holds your diaSources, their SNR,
and their parent obsHistId

2) the opsim table from the opsim run which generated your dia sources

fullerDiaSource.sql in the $MOPS_HACKS directory will create the
fullerDiaSource table.  If you're using Lynne's data with astrometric
error, I recommend you simply copy the contents from mops64's MySQL in
mops_noDeepAstromError.fullerDiaSource.  Or, you can simply run on
mops64.

None of the scripts will *write* to the database.



PREPARING THE TRACKLETS:

LinkTracklets uses bundled up sets of tracklets covering several days
of observation.  Our new script makeLinkTrackletsInput_byImages.py
will take care of creating scripts for running linkTracklets for you,
and can create input for C and/or C++ linkTracklets.

Before you run it, you'll need to make sure you have the needed tables
in your database, and organize your tracklets by start image, in order
to match the image obsHistIds to their observation times.

First:

$ mkdir byObsHistId/

For each tracklets file (e.g. 49524.maxv0.5.final.tracklets.byDiaIds)
you'll want to do this (use a shell script):

$ python $MOPS_HACKS binTrackletsByObsHist \
  49524.maxv0.5.final.tracklets.byDiaIds
  byObsHistId/

The script creates or appends to files like
byObsHistId/85679000.tracklets.byDiaId. 



MAKING THE LINKTRACKLETS INPUT FILES:

The script $MOPS_HACKS/makeLinkTracklets_byImages.py contains lots of
constants at the top, such as where to find the database, and where to
find the per-obsHist tracklets you created.  

There are also some additional parameters, like the size of the
tracking window for linkTracklets (15 is fast, 30 is slow).  You can
also set the max number of "start images" (images for which tracks are
started) per file (in case you want to break up the busy nights between
multiple runs).  

You also NEED to set the output directory for your linkTracklets
infiles: OUTPUT_LINKTRACKLETS_INFILE_DIR.

You can also set the flag WRITE_CPP_INFILES if you like to get input
for C++ linkTracklets.

Once you are content with the contents of the file, go ahead and run
the script with no arguments.  Then get ready for a healthy wait and a
lot of disk I/O; depending on the size of MAX_START_IMAGES_PER_FILE
this could take anywhere from 20 minutes to 10 hours or so.

When you're done, you shold have the following in OUTPUT_LINKTRACKLETS_INFILE_DIR:
For many Xs and Ys:

image_X_through_Y.info - human-readable info about the images in this data set, start images as well as potential second endpoint/support images.

image_X_through_Y.miti  - the tracklets, in MITI format.

image_X_through_Y.miti.diaIds - at each line N, the diaId of the
   detection at line N of the corresponding .miti file.

image_x_through_Y.start_t_range - the time offset of image Y - image X.



RUNNING LINKTRACKLETS (the new way, distributed by images)
--------------------------------------

Go to the directory where you created your .info, .miti,
.start_t_range files. Create a new directory called tracks/ and cd
there.  Then run:

$ python $MOPS_HACKS/makeLinkTrackletsRunScripts.py

This will create lots of .cmd.sh files in your current directory, one
for each of the .miti files in ../ .  Each one is a bash script which
runs linkTracklets for you with our curent set of favored parameters.

make sure to

$ setup auton

and you're ready to run linkTracklets - 

$ for CMD in *.cmd.sh 
  do
     bash $CMD
  done

Get ready to wait a few more hours!

This will create a lot of files ending in c.tracks.byIndices in the
current directory.  These are the tracks, in the
by-line-number-in-the-MITI-file-format.

You'll later want to convert these to byDiaId format.  Once you've
done so, you may find you have some redundant tracks, due to a known
quirk of C linkTracklets.  You can remove these via sort and uniq.

so for each image_X_through_Y.c.tracks.byIndices:

$ python $MOPS_HACKS/indicesToIds.py \
   image_X_through_Y.c.tracks.byIndices \
   ../image_X_through_Y.miti.diaIds
   image_X_through_Y.c.tracks.byDiaIds
$ cat image_X_through_Y.c.tracks.byDiaIds | sort -n | \
 uniq > image_X_through_Y.c.tracks.byDiaIds.unique




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

You'll need to tell the buildOrbitServerInput.py script a prefix for
these outfiles, e.g. myTracks.  It will create many files, e.g.

myTracks_0.in.request
myTracks_0.in.tracklet
myTracks_0.in.manifest
myTracks_1.in.request
... etc.

So in the directory where you want all these files written: 

$ python $MOPS_HACKS/buildOrbitServerInput.py /path/to/tracks.byDiaId myTracks



Next, to run orbit_server on all these files, use runOrbitServer.sh.
Modify the file so it points to the orbit_server.x executable on your
system, then just run it in the directory with the in.request,
in.tracklet, etc. files.


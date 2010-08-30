#/bin/csh

setlsst
setup mops_daymops
unsetup python

# this is the C++ findtracklets
setenv LINKTRACKLETS  $MOPS_DAYMOPS_DIR/bin/linkTracklets

# modify as needed 
setenv INPUT_DIR  $PWD/windows
setenv INPUT_BASE `ls $INPUT_DIR/*.dets | sed 's/.dets//g'`
setenv OUTPUT_DIR $PWD/tracks
if (! -e $OUTPUT_DIR) then
    mkdir $OUTPUT_DIR
endif

echo "Will run linkTracklets on : "
foreach INPUT ($INPUT_BASE)
    echo $INPUT
end
echo "And place output in directory " $OUTPUT_DIR

foreach INPUT ($INPUT_BASE)
    setenv BASE `basename $INPUT`
    setenv timelim1 `$BASE | sed 's/linkTracklets_input_//g'`
    set timelim = `echo $timelim1-0.5+.13 | bc`
    setenv CMD "$LINKTRACKLETS -d $INPUT.dets -t $INPUT.ids -o $OUTPUT_DIR/$BASE.tracks -L $timelim"
    echo "Running "$CMD
    #choose your method of timing  
    # gtime = GNU's version of time gives additional info if you have it - 
    #  it adds the memory usage, for example. 
    # available from http://directory.fsf.org/project/time/
    gtime --verbose --output=$OUTPUT_DIR/$BASE.time $CMD    
    # otherwise use system time command
    #time -o $OUTPUT_DIR/`basename $NIGHTLY .miti`.maxv0.5.time $CMD    
    echo ""
end

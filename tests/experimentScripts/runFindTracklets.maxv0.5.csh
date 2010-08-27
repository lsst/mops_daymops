#!/bin/csh

# this is the C++ findtracklets
setenv FINDTRACKLETS  $MOPS_DAYMOPS_DIR/bin/findTracklets

# modify as needed 
setenv NIGHTLY_DIASOURCES `ls $PWD/*.miti`
setenv OUTPUT_DIR $PWD/tracklets
if (! -e $OUTPUT_DIR) then
    mkdir $OUTPUT_DIR
endif

echo "Will run findTracklets on : "
foreach NIGHTLY ( $NIGHTLY_DIASOURCES)
    echo $NIGHTLY
end
echo "And place output in directory " $OUTPUT_DIR

foreach  NIGHTLY  ($NIGHTLY_DIASOURCES)
     setenv CMD "$FINDTRACKLETS -i $NIGHTLY -o $OUTPUT_DIR/`basename $NIGHTLY .miti`.maxv0.5.tracklets -v .5 -m 0.0"
    echo "Running "$CMD
    #choose your method of timing  
    # gtime = GNU's version of time gives additional info if you have it - 
    #  it adds the memory usage, for example. 
    # available from http://directory.fsf.org/project/time/
    gtime -o $OUTPUT_DIR/`basename $NIGHTLY .miti`.maxv0.5.time $CMD    
    # otherwise use system time command
    #time -o $OUTPUT_DIR/`basename $NIGHTLY .miti`.maxv0.5.time $CMD    
    echo ""
end

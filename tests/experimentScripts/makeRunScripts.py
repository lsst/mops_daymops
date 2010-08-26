#!/usr/bin/env python


"""

jmyers jun 10 2010

similar to the script for building linkTracklets PBS scripts, but this
time we make use of the list of start_t_range, so that we avoid
redundant processing of endpoint pairs which exist in multiple data sets.

Please be sure to set up auton first.

"""



import sys, os, os.path
from glob import glob


INPUT_TRACKLETS_DIR="/workspace1/jmyers/nightlyDiasNodeep/cppFindTracklets/15DayWindowsMaxv0.5/"
INPUT_TRACKLETS_SUFFIX=".tracklets.miti"

OUTPUT_SCRIPTS_DIR="/workspace1/jmyers/nightlyDiasNodeep/cppFindTracklets/15DayWindowsMaxv0.5/runScripts/"
SCRIPT_RESULTS_DIR="/workspace1/jmyers/nightlyDiasNodeep/cppFindTracklets/15DayWindowsMaxv0.5/runScripts/results/"

AUTON_DIR = os.getenv("AUTON_DIR")
LINKTRACKLETS=os.path.join(AUTON_DIR, "linkTracklets_modified/linkTracklets")
MAX_RA_ACC=.02
MAX_DEC_ACC=.02


MJD_FROM_LINE=lambda x : float(x.split()[1])



def firstTimeInFile(inputFile):
    minItem = None
    line = inputFile.readline()

    while line != "":
        lineMJD = MJD_FROM_LINE(line)        
        #print "minItem = ", minItem, ", lineMJD = ", lineMJD
        if minItem == None or minItem > lineMJD:
            minItem = lineMJD
            #print "MinItem was ", minItem, " changed to ", lineMJD
        line = inputFile.readline()

    return minItem


def getStartTimesFromInputFiles(inputTrackletsFiles):
    res = {}
    startTimes = []
    for f in inputTrackletsFiles:
        print "finding first tracklet start time in ", f, "..."
        startTime = firstTimeInFile(file(f,'r'))
        res[f] = startTime
        startTimes.append(startTime)
    return res, sorted(startTimes)



def getStartTRange(trackletFile, trackletFileToStartTimeMap, startTimes):
    """ given a current file, return the lowest start time of any file
    with a start time > the requested file.
    
    that is, if we have file a.tracklets which contains its earliest
    tracklet starting on 1234.56 and file b.tracklets with its
    earliest tracklet at 1244.13, then getStartTRange(a.tracklets)
    should return 1244.13.

    the idea being that we don't need to search for tracks starting on
    1244.13 or later when proccessing a.tracklets; those same tracks
    will be found when searching inside b.tracklets.
    """
    tStart = trackletFileToStartTimeMap[trackletFile]
    laterTimes = filter(lambda x: x > tStart, startTimes) 
    if laterTimes == []:
        return None
    else:
        return min(laterTimes)
    



def mkRunscript(trackletFile, trackletFileToStartTimeMap, startTimes):
    outputBN=SCRIPT_RESULTS_DIR + os.path.basename(trackletFile)[:-(len(INPUT_TRACKLETS_SUFFIX))]
    startTRange =  getStartTRange(trackletFile, trackletFileToStartTimeMap, startTimes) 

    if startTRange == None:
        startTRangeStr = ""
    else:
        startTRangeStr = " start_t_range " + str(startTRange)

    s = "#!/bin/bash" 
    s += "\n"
    s += "\n"
    s += "CMD=\"" + LINKTRACKLETS + " file " + trackletFile + " indicesfile " + outputBN + ".c.tracks.byIndices " +\
    " acc_r " + str(MAX_RA_ACC) + " acc_d " + str(MAX_DEC_ACC) + startTRangeStr + "\"\n"
    s += """echo Running $CMD...
date
/usr/bin/time -o """ + outputBN + """.time $CMD | tee """ + outputBN + """.log
echo Finished at `date`
"""
    return s


if __name__=="__main__":
    import pickle

    inputTracklets = glob(INPUT_TRACKLETS_DIR + "*" + INPUT_TRACKLETS_SUFFIX)
    print inputTracklets

    try:
        trackletFileToStartTimeMap = pickle.load(file("trackletFileToStartTimeMap.pickle", 'r'))
        allStartTimes = pickle.load(file("allStartTimes.pickle",'r'))
    
        print "Loaded start time information from pickle files."

    except:
        print "Calculating start time information..."
        trackletFileToStartTimeMap, allStartTimes = getStartTimesFromInputFiles(inputTracklets)
    
        pickle.dump(trackletFileToStartTimeMap, file("trackletFileToStartTimeMap.pickle",'w'))
        pickle.dump(allStartTimes, file("allStartTimes.pickle",'w'))


    print trackletFileToStartTimeMap
    print allStartTimes

    for trackletFile in trackletFileToStartTimeMap:
        script = mkRunscript(trackletFile, trackletFileToStartTimeMap, allStartTimes)
        print script
        outfile = file(OUTPUT_SCRIPTS_DIR+os.path.basename(trackletFile[:-(len(INPUT_TRACKLETS_SUFFIX))]),'w')
        outfile.write(script)
                       

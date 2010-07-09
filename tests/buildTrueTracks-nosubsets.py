#!/usr/bin/env python


"""

jmyers 

2/11/10

read through detections with ground-truth data.  We ASSUME that the
data is sorted first by object ID.

Find the true tracks that would be found by
findTracklets/linkTracklets and write them to output.

Caveats:

Cheats somewhat on guessing night numbers, does not use GMT. 

When simulating linkTracklets, does not make use of the quality-of-fit
of detections to endpoints: only the acceleration implied by the
endpoint tracklets is considered when deciding if tracklets should be
used to form a track.

FindTracklets simulation does not attempt to simulate the correct
behavior of findTracklets (PHT mode) or
findTracklets+CollapseTracklets. Instead, we simply generate at most
one tracklet per detection.

"""

import numpy
import time
import os

DEBUG=False



class Detection(object):
    def __init__(self, time, ra, dec):
        self.time = time
        self.ra = ra
        self.dec = dec



# TBD: Need to replace this with something that uses GMT...
def nightNumFromMJD(mjd):
    guess = int(mjd)
    if (mjd - guess) < .5:
        return guess
    else:
        return guess + 1




def isFindableFromTimes(times, minPerNight=2):
    """                                                                                                                                                                                                                          
    returns true if on at least three distinct nights, there are >=                                                                                                                                                              
    minPerNight items in 'times'                                                                                                                                                                                                 
    """
    nightHistogram = {}
    for time in times:
        nn = nightNumFromMJD(time)
        if nightHistogram.has_key(nn):
            nightHistogram[nn] += 1
        else:
            nightHistogram[nn] = 1

    usableNights = []
    for night in nightHistogram.keys():
        if nightHistogram[night] >= minPerNight:
            usableNights.append(night)

    if len(usableNights) >= 3:
        return True
    else:
        return False




def writeTracks(objectDets, velocityLimit, trackletTimeLimit,
                trackletMinTime,
                raAccelerationLimit, decAccelerationLimit,
                linkTrackletsTimeWindow, outfile):
    """
    objectDets: a list of Detections attributable to the SAME, SINGLE
    object.

    velocityLimit is the limit on velocity for two detections to be
    considered as part of the same tracklet.

    trackletTimeLimit is the limit on time difference between
    exposures for two detections to be linked into a tracklet.

    ra/dec Acceleration limit are the maximum amount of acceleration
    between two endpoint tracklets for a tracklet to be formed.

    linkTrackletsTimeWindow is the maximum time between two tracklets
    for them to be considered in a track.


    ------

    find all "linkable" detections under these constraints and write
    them to a file.  Require at least 6 detections, with tracklets on
    at least 3 nights (detections are considered as from separate
    nights if they are >= .5 days apart).  Write these to outfile.

    return the number of tracks generated for this object.

    apr. 27: also, return tracklet start times for all tracklets in
    some track
    """
    
    if DEBUG:
        print "Saw another object with detections:"
        for i in objectDets:
            print i.time, i.ra, i.dec
        print ""

    tracklets = findTracklets(objectDets, velocityLimit, trackletTimeLimit, trackletMinTime)

    if DEBUG:
        if len(tracklets) == 0:
            print "Got no tracklets!"
        for tracklet in tracklets:
            p1 = objectDets[tracklet[0]]
            p2 = objectDets[tracklet[1]]
            print "\tGot tracklet: "
            print "\t",p1.time, p1.ra, p1.dec
            print "\t",p2.time, p2.ra, p2.dec
            print ""
    
    tracks = linkTracklets(objectDets, tracklets, raAccelerationLimit, decAccelerationLimit,
                           linkTrackletsTimeWindow)
    if DEBUG: 
        print "got tracks: ", tracks
        if (len(tracks) == 0):
            print "Could not find object with dets: "
            for det in objectDets:
                print det.time, det.ra, det.dec

    #TBD: write to output...?

    trackletStartTimes = []
    for track in tracks:
        thisTrackStartTimes = []
        for tracklet in track:
            thisTrackStartTimes.append(objectDets[tracklet[0]].time)
        trackletStartTimes.append(thisTrackStartTimes)

    #if len(tracks) > 0 :
    #    print tracks
    #    print trackletStartTimes
    return len(tracks), trackletStartTimes



def trackletRaVelocity(objectDets, tracklet):
    return (objectDets[tracklet[0]].ra - objectDets[tracklet[1]].ra) / \
        (objectDets[tracklet[0]].time - objectDets[tracklet[1]].time)




def trackletDecVelocity(objectDets, tracklet):
    return (objectDets[tracklet[0]].dec - objectDets[tracklet[1]].dec) / \
        (objectDets[tracklet[0]].time - objectDets[tracklet[1]].time)





def linkTracklets(objectDets, tracklets, raAccLimit, decAccLimit, timeWindow):
    """
    Using detections which are known to belong to the same object,
    give back the same tracks that linkTracklets would generate (or as
    close an approximation as we can get)
    """
    tracks = []

    #consider every tracklet as a starting point, just like real linkTracklets
    for i in range(len(tracklets)):

        if DEBUG == True:
            print "Looking for tracks starting at tracklet #", i, " of ", len(tracklets)

        startTracklet = tracklets[i]
        startTime = objectDets[startTracklet[0]].time

        startRaV = trackletRaVelocity(objectDets, startTracklet)
        startDecV = trackletDecVelocity(objectDets, startTracklet)

        timeCompatibleTracklets = [ ]

        #print "all tracklets: "
        #print tracklets

        #find all tracklets which are compatible in the time domain.
        
        for j in range(i + 1, len(tracklets)):
            secondEndpointTime = objectDets[tracklets[j][0]].time
            if secondEndpointTime - startTime <= timeWindow and\
                    secondEndpointTime - startTime > 0.:
                # we will later filter on whether there were detections on
                # >=3 nights; for now only worry about the time window.
                timeCompatibleTracklets.append(tracklets[j])

        #if DEBUG:
        #    print "time-compatible tracklets:"
        #    for t in timeCompatibleTracklets:
        #        print [[objectDets[t[z]].time, objectDets[t[z]].ra, objectDets[t[z]].dec] for z in [0,1]]
                
        accCompatibleTracklets  = [ ]

        # among those which are compatible in time, check if they are compatible
        # in acceleration.
        for j in range(len(timeCompatibleTracklets)):
            jRaV = trackletRaVelocity(objectDets, timeCompatibleTracklets[j])
            jDecV = trackletDecVelocity(objectDets, timeCompatibleTracklets[j])
            jTime = objectDets[timeCompatibleTracklets[j][0]].time

            if abs(jRaV - startRaV) / abs(jTime - startTime) <= raAccLimit:
                impliedAcc = abs(jDecV - startDecV) / abs(jTime - startTime)
                if impliedAcc <= decAccLimit:
                    accCompatibleTracklets.append(timeCompatibleTracklets[j])
                #else:
                #    print "found object accelerating too quickly: %f deg/day/day" % impliedAcc
        #if DEBUG:
        #    print "time + acc compatible tracklets:"
        #    for t in accCompatibleTracklets:
        #        print [[objectDets[t[z]].time, objectDets[t[z]].ra, objectDets[t[z]].dec] for z in [0,1]]

        #new version: 
        # only find maximal tracks. Since we don't bother to account for 
        # quality of fit to the best-fit quadratic, we can cheat a bit and do this quickly.
        # This will be much faster.
        trackletStartTimes = []
        trackSubset = [startTracklet]
        for j in range(len(accCompatibleTracklets)):            

            jTime = objectDets[accCompatibleTracklets[j][0]].time
            trackSubset.append(accCompatibleTracklets[j])
            trackletStartTimes.append(jTime)

        if DEBUG:
            print "Built maximal track, got tracklet start times: ", trackletStartTimes                

        #isFindableFromTimes checks whether we have detections from three distinct nights.
        # if we do, add them to output.
        if isFindableFromTimes(trackletStartTimes, minPerNight=1):
            if DEBUG:
                print "Track spans at least three nights, adding to results."
            tracks.append(trackSubset)
        else:
            if DEBUG:
                print "Track does not span sufficient nights. abandoning it."

        #print "Current set of findable tracks starting at", tracklets[i], ":"
        #for t in tracks:
        #    print t

    #print "Final set of tracks for this object:"
    #for t in tracks:
    #    print t
    return tracks

            


def euclideanDistance(detection1, detection2):
    dRa = detection1.ra - detection2.ra
    dDec = detection1.dec - detection2.dec
    return numpy.sqrt(dRa**2 + dDec**2)





def findTracklets(objectDets, velocityLimit, trackletTimeLimit, trackletMinTime):
    output = []
    for i in range(len(objectDets)):
        #shortcut to prevent n^2 growth of tracklets in "deep stacks":
        # allow only one tracklet per detection. 
        haveTracklet = False
        for j in range(i + 1, len(objectDets)):
            if not haveTracklet:
                dTime = abs(objectDets[i].time - objectDets[j].time)
                if dTime <= trackletTimeLimit and dTime >= trackletMinTime:
                    trackletVelocity = euclideanDistance(objectDets[i], objectDets[j])/dTime
                    if trackletVelocity <= velocityLimit:
                        output.append([i,j])
                        haveTracklet = True
    return output



def increment(d, key):
    if d.has_key(key):
        d[key] += 1
    else:
        d[key] = 1


def inRangeNonInclusive(t, a, b):
    return t > a and t < b


if __name__ == "__main__":
    import sys

    [velocityLimit, trackletTimeLimit, trackletMinTime] = map(float, sys.argv[1:4])
    [raAccelerationLimit, decAccelerationLimit] = map(float, sys.argv[4:6])
    linkTrackletsTimeWindow = float(sys.argv[6])
    [infile, outfile] = sys.argv[7:]

    infile = open(infile, 'r')
    outfile = open(outfile, 'w')

    print "Got tracklet velocity limit:              %f" % velocityLimit
    print "Got tracklet min time lmit:               %f" % trackletMinTime
    print "Got tracklet time limit:                  %f" % trackletTimeLimit
    print "Got track Ra Acc. limit:                  %f" % raAccelerationLimit
    print "Got trac Dec Acc. limit:                  %f" % decAccelerationLimit
    print "Got track time window:                    %f" % linkTrackletsTimeWindow

    line = infile.readline().split()

    curObjectData = [ ]

    nObjects = 0
    nTrackGeneratingObjects = 0

    imgTimesToTrackCounts = {}
    imgTimesToLastEndpointCounts = {}
    imgTimesToFirstEndpointCounts = {}
    
    #timeWindowDicts[x][y] == number of tracks which start on night y and have a valid second endpoint on night x
    timeWindowDicts = {10:{}, 11:{}, 12:{}, 13:{}, 14:{}, 15:{}, 16:{}, 17:{}, 18:{}, 19:{}, 20:{}}

    statusLog = file('/tmp/buildTrueTracs-nosubsets.py' + str(os.getpid()) + '_progress.log', 'a')
    while line != []:

        curObject = line[0]
        [mjd, ra, dec] = map(float, line[1:4])
        d = Detection(mjd, ra, dec)
        curObjectData.append(d)

        nextLine = infile.readline().split()
        done = False

        #check whether this was the last line containing this object's data
        if nextLine == []:
            done = True
        elif nextLine[0] != curObject:
            done = True

        if done:
            #we've finished reading data for one object into curObjectData.
            nObjects += 1
            statusLog.write("working on object " + str(nObjects) + "\n")
            statusLog.flush()
            print "Working on object", nObjects
            LOCALDEBUG=False
            if nObjects == 1855:
                LOCALDEBUG=True

            nTracks, allTrackletStartTimes = writeTracks(curObjectData, velocityLimit, trackletTimeLimit,
                                                         trackletMinTime,
                                                         raAccelerationLimit, decAccelerationLimit,
                                                         linkTrackletsTimeWindow, outfile)
            if LOCALDEBUG:
                print "Got ", len(allTrackletStartTimes), " tracks."
            for trackletStartTimes in allTrackletStartTimes:                
                if LOCALDEBUG:
                    print "Working on track with ", len(trackletStartTimes), " tracklets."
                trackletsSoFar = 0
                for trackletStartTime in trackletStartTimes:                
                    trackletsSoFar += 1
                    increment(imgTimesToTrackCounts, trackletStartTime)
                    nTracklets = len(trackletStartTimes)

                    if nTracklets < 3:
                        print "SOMETHING WEIRD: LESS THAN THREE TRACKLETS IN TRACK!"
                    else:
                        if trackletsSoFar >= 2:
                            increment(imgTimesToLastEndpointCounts, trackletStartTime)
                        if trackletsSoFar <= nTracklets - 2:
                            increment(imgTimesToFirstEndpointCounts, trackletStartTime)
                        
                trackletStartNights = sorted(map(nightNumFromMJD,  trackletStartTimes))
                for night in trackletStartNights:                                                
                    for offset in range(10,21):
                        #print "night + offset = ", night+offset
                        if night + offset in trackletStartNights:
                            #print "     found track starting at", night, " ending at ", night + offset
                            intermediateNights = filter(lambda x: inRangeNonInclusive(x, night, night + offset), trackletStartNights)
                            #print "      intermediate nights were ", intermediateNights
                            if (len(intermediateNights)) > 0:
                                #print "      it DID have an intermediate night"
                                increment(timeWindowDicts[offset], night)
                        
            curObjectData = []
            if nTracks > 0:
                nTrackGeneratingObjects += 1
            #if nObjects % 1 == 0:
            #    print time.asctime()
            #    print "SO FAR, of %i objects, %i generated tracks (%f%%)" % \
            #        (nObjects, nTrackGeneratingObjects, nTrackGeneratingObjects*100./nObjects)


        #get the next line        
        line = nextLine

        
        
    print "of %i objects, %i generated tracks (%f%%)" % \
        (nObjects, nTrackGeneratingObjects, nTrackGeneratingObjects*100./nObjects)
        
    print "Image times with viable start tracklets for some track: "
    for imgTime in sorted(imgTimesToTrackCounts.keys()):
        print imgTime, " : ", imgTimesToTrackCounts[imgTime], " tracks contained a tracklet rooted in this image."

        if imgTimesToFirstEndpointCounts.has_key(imgTime):
            print "\t", imgTimesToFirstEndpointCounts.has_key(imgTime), " tracks could start with a tracklet rooted in this image."

        if imgTimesToLastEndpointCounts.has_key(imgTime):
            print "\t", imgTimesToLastEndpointCounts.has_key(imgTime), " tracks could end with a tracklet rooted in this image."

    print ""
    print ""

    for timeWindow in range(10, 20):
        for startNight in sorted(timeWindowDicts[timeWindow].keys()):
            print "There are ", timeWindowDicts[timeWindow][startNight], " plausible tracks starting on night ", startNight, " and ending on night ", startNight + timeWindow
        print ""

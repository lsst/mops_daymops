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

DEBUG=False
OUTPUT_FORMAT="diaIds" # "miti" or "diaIds"
OBSCODE="807"
FORCED_MAG=24.0
CUR_TRACK_ID=0 #this is the first track ID to be written, successive tracks will have successive IDs...



class Detection(object):
    def __init__(self, time, ra, dec, detId=None, objId=None):
        self.time = time
        self.ra = ra
        self.dec = dec
        self.detId = detId
        self.objId = objId



# TBD: Need to replace this with something that uses GMT...
def nightNumFromMJD(mjd):
    return int(mjd)





def isFindableFromTimes(times, minPerNight=2, minNights=3):
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

    if len(usableNights) >= minNights:
        return True
    else:
        return False




def writeTracks(objectDets, velocityLimit, trackletTimeLimit,
                trackletMinTime,
                raAccelerationLimit, decAccelerationLimit,
                linkTrackletsTimeWindow, outfile, findableOutfile, curTrackId,
                minObs, minNights):
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
                           linkTrackletsTimeWindow, minObs, minNights)
    if DEBUG:
        print "got tracks: ", tracks
        if (len(tracks) == 0):
            print "Could not find object with dets: "
            for det in objectDets:
                print det.time, det.ra, det.dec

    #write tracks for this object...
    curTrackId = writeTracksToFile(objectDets, tracks, tracksOutfile, curTrackId)

    if len(tracks) > 0:
        objName = objectDets[tracks[0][0][0]].objId
        findableOutfile.write(objName + "\n")

    trackletStartTimes = []
    for track in tracks:
        for tracklet in track:
            trackletStartTimes.append(objectDets[tracklet[0]].time)
    return len(tracks), trackletStartTimes, curTrackId



def trackletRaVelocity(objectDets, tracklet):
    return (objectDets[tracklet[0]].ra - objectDets[tracklet[1]].ra) / \
        (objectDets[tracklet[0]].time - objectDets[tracklet[1]].time)




def trackletDecVelocity(objectDets, tracklet):
    return (objectDets[tracklet[0]].dec - objectDets[tracklet[1]].dec) / \
        (objectDets[tracklet[0]].time - objectDets[tracklet[1]].time)




def writeTracksToFile(objDets, tracks, tracksOutfile, trackId):
    for track in tracks:
        detIds = []
        objName = None
        for tracklet in track:
            for detIndex in tracklet:
                detIds.append(objDets[detIndex].detId)
                newObjName = objDets[detIndex].objId
                if objName is None:
                    objName = newObjName
                if objName != newObjName:
                    print "ERROR: in writeTracksToFile: found invalid track?!"
                    print "Cur detection is from ", newObjName, " other detection was from ", objName
                    print "All dets so far are: ", detIds
                    sys.exit(1)
        if OUTPUT_FORMAT=="diaIds":
            for detId in detIds:
                tracksOutfile.write("%s " % detId)
            tracksOutfile.write('\n')
            #tracksToObjNameFile.write("%s\n" % objName)
        elif OUTPUT_FORMAT=="miti":
            for tracklet in track:
                for detIndex in tracklet:
                    det = objDets[detIndex]
                    # ID EPOCH_MJD RA_DEG DEC_DEG MAG OBSCODE OBJECT_NAME LENGTH ANGLE [EXPOSURE_TIME]
                    tracksOutfile.write("%d %f %f %f %f %s %s 0.0 0.0\n" % \
                                            (trackId, det.time, det.ra, det.dec, FORCED_MAG, 
                                             OBSCODE, det.objId))
        else:
            raise Exception("unknown output format %s!" % OUTPUT_FORMAT)

    trackId += 1
    return trackId


def linkTracklets(objectDets, tracklets, raAccLimit, decAccLimit, timeWindow, minObs, minNights):
    """
    Using detections which are known to belong to the same object,
    give back the same tracks that linkTracklets would generate (or as
    close an approximation as we can get)
    """
    tracks = []

    #consider every tracklet as a starting point, just like real linkTracklets
    for i in range(len(tracklets)):

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

        if DEBUG:
            print "time-compatible tracklets:"
            for t in timeCompatibleTracklets:
                print [[objectDets[t[z]].time, objectDets[t[z]].ra, objectDets[t[z]].dec] for z in [0,1]]
                
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
        if DEBUG:
            print "time + acc compatible tracklets:"
            for t in accCompatibleTracklets:
                print [[objectDets[t[z]].time, objectDets[t[z]].ra, objectDets[t[z]].dec] for z in [0,1]]

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
            print "Built maximal track, got tracklet start times: ", [startTime] + trackletStartTimes

        #isFindableFromTimes checks whether we have detections from three distinct nights.
        # if we do, add them to output.
        numObs = 0
        for tracklet in trackSubset:
            for dia in tracklet:
                numObs += 1
        if numObs >= minObs:
            
            if isFindableFromTimes(trackletStartTimes, minPerNight=1, minNights=minNights):
                tracks.append(trackSubset)
            else:
                if DEBUG:
                    print "Discarded a track which is not observed on sufficient number of nights."
        else:
            if DEBUG:
                print "Discarded a track with ", numObs, " observations < ", minObs

        #print "Current set of findable tracks starting at", tracklets[i], ":"
        #for t in tracks:
        #    print t

    
    return tracks

            


def euclideanDistance(detection1, detection2):
    dRa = detection1.ra - detection2.ra
    dDec = detection1.dec - detection2.dec
    return numpy.sqrt(dRa**2 + dDec**2)





def findTracklets(objectDets, velocityLimit, trackletTimeLimit, trackletMinTime):
    output = []
    isInTrackletAlready = []
    for i in range(len(objectDets)):
        isInTrackletAlready.append(False)
    for i in range(len(objectDets)):
        # simulate "deep stacks" or whatever by clumping together all
        # detections from the same night into one giant tracklet

        #this isn't quite ideal but it should be good enough, I guess...
        if not isInTrackletAlready[i]:
            if DEBUG:
                print "Attempting to build a tracklet starting at ", objectDets[i].time, objectDets[i].ra, objectDets[i].dec
            newTracklet = [i]
            isInTrackletAlready[i] = True
            for j in range(i + 1, len(objectDets)):
                dTime = abs(objectDets[i].time - objectDets[j].time)
                if DEBUG:
                    print "Attempting to add detection", objectDets[j].time, objectDets[i].ra, objectDets[i].dec, ", dtime = ", dTime, " so far tracklet is size ", len(newTracklet)
                    
                if (dTime <= trackletTimeLimit and dTime >= trackletMinTime) or (len(newTracklet) >= 2 and dTime <= .5):
                    trackletVelocity = euclideanDistance(objectDets[i], objectDets[j])/dTime
                    if DEBUG:
                        print "Velocity was ", trackletVelocity
                    if trackletVelocity <= velocityLimit:
                        if DEBUG:
                            print "Added it to tracklet."
                        output.append([i,j])
                        isInTrackletAlready[j] = True
    return output




if __name__ == "__main__":
    import sys

    [velocityLimit, trackletTimeLimit, trackletMinTime] = map(float, sys.argv[1:4])
    [raAccelerationLimit, decAccelerationLimit] = map(float, sys.argv[4:6])
    minObs, minNights = map(int, sys.argv[6:8])
    linkTrackletsTimeWindow = float(sys.argv[8])
    [infile, tracksOutfile, findableObjects] = sys.argv[9:]

    infile = open(infile, 'r')
    tracksOutfile = open(tracksOutfile, 'w')
    findableOutfile = open(findableObjects, 'w')

    print "Got tracklet velocity limit:              %f" % velocityLimit
    print "Got tracklet min time lmit:               %f" % trackletMinTime
    print "Got tracklet time limit:                  %f" % trackletTimeLimit
    print "Got track Ra Acc. limit:                  %f" % raAccelerationLimit
    print "Got trac Dec Acc. limit:                  %f" % decAccelerationLimit
    print "Got track time window:                    %f" % linkTrackletsTimeWindow
    
    print "Writing tracks, as sets of det IDs, to    %s" % tracksOutfile
    print "Writing findable object IDs to            %s" % findableOutfile

    line = infile.readline().split()

    curObjectData = [ ]

    nObjects = 0
    nTrackGeneratingObjects = 0
    curTrackId = 0

    imgTimesToTrackCounts = {}

    while line != []:


        #modified version uses MITI input
        curObject = line[6]
        mjd = line[1]
        ra = line[2]
        dec = line[3]
        detId = line[0]
        # old format: uses output of 
        #select diaSourceId, taiMidPoint, ra, decl, groundTruthMovingObjectId
        #curObject = line[4]
        #mjd = float(line[1])
        #ra = float(line[2])
        #dec = float(line[3])
        #detId = int(line[0])
        
        d = Detection(float(mjd), float(ra), float(dec), detId=detId, objId=curObject)
        curObjectData.append(d)

        nextLine = infile.readline().split()
        done = False

        #check whether this was the last line containing this object's data
        if nextLine == []:
            done = True
        #elif nextLine[4] != curObject:
        elif nextLine[6] != curObject:
            done = True

        if done:
            #we've finished reading data for one object into curObjectData.
            nObjects += 1
            #if nObjects == 48:
            #    DEBUG=True
            #print "Found object with ", len(curObjectData), " detections."
            #for data in curObjectData:
            #    print "Got object with ID, MJD, ra, dec, detId: "
            #    print data.objId, data.time, data.ra, data.dec, data.detId
            #sys.exit(1)
            

            nTracks, trackletStartTimes, curTrackId = writeTracks(curObjectData, velocityLimit, trackletTimeLimit,
                                                                  trackletMinTime,
                                                                  raAccelerationLimit, decAccelerationLimit,
                                                                  linkTrackletsTimeWindow, tracksOutfile, findableOutfile,
                                                                  curTrackId, minObs, minNights)
            for trackletStartTime in trackletStartTimes:
                if imgTimesToTrackCounts.has_key(trackletStartTime):
                    imgTimesToTrackCounts[trackletStartTime] += 1
                else:
                    imgTimesToTrackCounts[trackletStartTime] = 1

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
        print imgTime, " : ", imgTimesToTrackCounts[imgTime]

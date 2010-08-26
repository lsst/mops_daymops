#!/usr/bin/env python

"""

jmyers jul. 6 2010

this is a script for testing the ranging routines provided by
Francesco.  Provide some tracks (sets of detections) and it will write out the
resulting ranging orbits.

"""


DEFAULT_OBSCODE='807'

import numpy, oorb

class Detection(object):
    def __init__(self, ra, dec, mjd, mag, objId):
        self.ra = ra
        self.dec = dec
        self.mjd = mjd
        self.mag = mag
        self.objId = objId

    def __init__(self, mitiString):
        items = mitiString.strip().split()
        self.mjd, self.ra, self.dec, self.mag = map(float, items[1:5])
        self.objId = items[6]



class Track(object):
    def __init__(self, detections):
        self.detections = detections
        self.objectName = None
    def addDetection(self, det):
        self.detections.append(det)
    def getObjectName(self):
        return self.objectName
    def setObjectName(self, name):
        self.objectName = name

    def isTrue(self):
        objName = None
        isTrue = True
        for det in self.detections:
            if objName == None:
                objName = det.objId
            if det.objId != objName or det.objId == 'NS':
                return False
        return True


class Orbit(object):
    def __init__(self, q, e, i, node, argPeri, meanAnom, epoch, src):
        self.q=q
        self.e=e
        self.i=i
        self.node=node
        self.argPeri=argPeri
        self.meanAnom=meanAnom
        self.epoch=epoch
        self.src=src

def orbitDetermination(track, 
                       elementType='keplerian', 
                       numRangingOrbits=5000,
                       stdDev=8.3333333333333331e-05,
                       obscode='566'):
    """
    jmyers - this is stolen from Francesco's linking.py in r.13810.  I
    had to hack on it a bit because we're not tethered to a DB and
    this code assumes that tracks are sets of tracklets, which isn't
    strictly true.  It also expected fluxes instead of mags; we have
    mags in MITI format.
    """
    # We do statistical ranging first. Then we do LSL on the orbits we get from 
    # statistical ranging.
    trackId = 0
    coords = []
    mjds = []
    mags = []
    filters = []
    obscodes = []                           # In reality we only have one!
    
    # Get to the DiaSources.
    for d in track.detections:
        coords.append([d.ra, stdDev, d.dec, stdDev])
        mags.append(d.mag)
        mjds.append(d.mjd)
        # TODO: do we need the filter name or is the ID string OK?
        filters.append('0') #jmyers - hacking this, but does it matter?
            
    # Now convert those to numpy arrays.
    coords = numpy.array(coords, dtype='d')
    coords.shape = (len(coords), 4)
    obscodes = [obscode, ] * len(mjds)
    mjds = numpy.array(mjds, dtype='d')
    mags = numpy.array(mags, dtype='d')

    
    # We can start statistical ranging.
    try:
        # Choose just 3 or 4 detections for ranging.
        rangingOrbits = oorb.ranging_fast(trackId=0,
                                          coords=coords[:3],
                                          mjds=mjds[:3],
                                          mags=mags[:3],
                                          obscodes=obscodes[:3],
                                          filters=filters[:3],
                                          elementType=elementType,
                                          numOrbits=numRangingOrbits)

        # Now pass them to LSL and this time use all detections.
        # res = (out_orbit, out_covariance, out_sigmas, out_correlation)
        res = oorb.lsl_fast(trackId=trackId,
                            coords=coords,
                            mjds=mjds,
                            mags=mags,
                            obscodes=obscodes,
                            filters=filters,
                            rangingOrbits=rangingOrbits)
    except:
        #Orbit determination failed for this track. Oh well.
        return(None)
    
    # If everything went well, we have an orbit with covariance.
    # res[0]: [trackId, a, e, i, node, argPeri, m, epoch, H, G, elTypeId]
    # res[1]: is a 6x6 covariance matrix, get the upper diagonal form.
    cov = []
    for i in (0, 1, 2, 3, 4, 5):
        for j in range(i, 6, 1):
            cov.append(res[1][i][j])
    
    # q = (1 - e) * a
    [trackId, 
     a, 
     e, 
     i, 
     node, 
     argPeri, 
     m, 
     epoch, 
     H, 
     G, 
     elTypeId] = list(res[0])
    return(Orbit(q=(1. - e) * a,
                 e=e,
                 i=i,
                 node=node,
                 argPeri=argPeri,
                 meanAnom=m,
                 epoch=epoch,
                 src=cov))




def getTracksFromMitiFile(lines):
    """ given a set of tracks in MITI format (that is, a set of
    detections where detections from the same track share the same ID)
    return a list of tracks (which themselves are lists of
    detections)"""

    tracks = {}

    for line in lines:
        items = line.strip().split()
        if len(items) < 9 or len(items) > 10:
            raise Exception("are input data file lines in MITI format?  Found one had %i items" % len(items))
        curTrackId = int(items[0])

        if not tracks.has_key(curTrackId):
            tracks[curTrackId] = Track([])

        newDet = Detection(line)
        tracks[curTrackId].addDetection(newDet)
        

    return tracks




def getTracksFromDetsAndIndicesFiles(detsInFile, indicesInfile):
    tracks = {}
    curTrackId = 0
    dets = detsInFile.readlines()
    for line in indicesInfile.readlines():
        indices = map(int, line.strip().split())

        newTrack = Track([])

        for index in indices:
            det = Detection(dets[index])
            newTrack.addDetection(det)

        tracks[curTrackId] = newTrack
        curTrackId += 1
    return tracks





def writeDeltaToFile(dbc, deltasOutfile, orbit, objId, trackId):
    if (objId != 'UNTRUE_TRACK'):
        dbc.execute(""" select q,e,i,node, argperi, epoch from orbits where ssmId=%d; """ % (int(objId)))
        res = dbc.fetchone()
        [q, e, i, node, argPeri, epoch] = res
        print>>deltasOutfile, "objId = ", objId, "trackId = ", trackId, epoch, \
            q - orbit.q, e - orbit.e, i - orbit.i, node - orbit.node, argPeri - orbit.argPeri
    else:
        print>>deltasOutfile, "UNTRUE TRACK ", trackId, " got a successful orbit!"



def writeFailureToDeltaFile(deltasOutfile, objId, trackId):
    print>>deltasOutfile, "objId =", objId, "trackId =", trackId, " IOD FAILED!"




if __name__=="__main__":
    import sys,time
    import MySQLdb

    db = MySQLdb.connect(user="jmyers", passwd="jmyers", db="newOrbitsFeb2010")
    dbc = db.cursor()

    if len(sys.argv) == 4:
        print "Reading tracks in MITI format (ID of each detection is taken to be ID of parent track) from %s" % sys.argv[1]
        infile = file(sys.argv[1],'r')
        outfile = file(sys.argv[2],'w')
        deltasOutfile = file(sys.argv[3], 'w')
        tracks = getTracksFromMitiFile(infile.readlines())

    if len(sys.argv) == 5:
        print "Reading tracks from sets of MITI detections (%s) and a set of sets of indices into that file (%s) " % (sys.argv[1], sys.argv[2])
        detsInfile = file(sys.argv[1], 'r')
        indicesInfile = file(sys.argv[2], 'r')
        outfile = file(sys.argv[3],'w')
        deltasOutfile = file(sys.argv[4], 'w')
        tracks = getTracksFromDetsAndIndicesFiles(detsInfile, indicesInfile)
    else:
        print """takes either one input argument (MITI file of tracks)
        or two input arguments (MITI file of detections, file of
        detection indices which form the tracks) and two output
        parameters: one for orbits, one for deltas"""

        sys.exit(1)

    count = 0
    for id in tracks.keys():
        if tracks[id].isTrue():
            print "Track ", id, " : "
            dets = tracks[id].detections
            for det in dets:
                print "\t", det.mjd, det.ra, det.dec, det.mag, det.objId

    rangingTime = 0
    rangingCalls = 0
    rangingSuccess = 0
    print "Calling ranging..."
    for trackId in tracks.keys():
        #if tracks[trackId].isTrue():
        print "\trunning on track ", trackId
        print "\ttime is ", time.asctime()
        t0 = time.time()
        orbit = orbitDetermination(tracks[trackId], obscode=DEFAULT_OBSCODE, stdDev=.001)
        rangingTime = time.time() - t0
        rangingCalls += 1

        if tracks[trackId].isTrue():
            objId = tracks[trackId].detections[0].objId
        else:
            objId = 'UNTRUE_TRACK'

        if orbit != None:
            rangingSuccess += 1
            print>>outfile, "Track ID = ", trackId, " corresponds to object: ", objId, "; fit orbit: ", orbit.q, orbit.e, orbit.i, orbit.node, orbit.argPeri, orbit.meanAnom, orbit.epoch, orbit.src
            print "\tDone! Got an orbit."
            print "Looking up orbit in DB and writing deltas..."                
            writeDeltaToFile(dbc, deltasOutfile, orbit, objId, trackId)
            print "Done!"
        else:
            print "\tOrbit determination FAILED for track ", trackId
            print "\tTime is ", time.asctime()
            print>>outfile, "Track ID = ", trackId, " corresponds to object: ", objId, " FAILED ORBIT DETERMINATION"
            writeFailureToDeltaFile(deltasOutfile, objId, trackId)
        print "\ttime is ", time.asctime(), "\n"
        count += 1
        if (count > 100):
            print "Spent %f sec on %d calls to ranging. %d succeeded." % (rangingTime, rangingCalls, rangingSuccess)

    print "Spent %f sec on %d calls to ranging. %d succeeded." % (rangingTime, rangingCalls, rangingSuccess)

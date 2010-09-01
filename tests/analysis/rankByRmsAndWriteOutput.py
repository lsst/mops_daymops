#!/usr/bin/env python

""" jmyers aug 30 2010


Takes the infile and outfile from C linkTracklets (i.e. tracklets in
MITI format, a set of indices into that file) and then finds the N
"best" tracks associated with each detection.  Output some stats about
how many true tracks survived and how many total tracks were saved,
etc.

"best"-ness will be determined like this:

1) a track which is a superset of another is always better (replace
the shorter one) 
2) after that, lower RMS -> better.


As this script got extended, it now also writes the following for each track it sees:
 - track line number from the tracks file
 - whether it is a true track
 - length of the track
 - RMS of the track
 - the track's initial declination
 - the fit error for EACH detection

"""

TRACKS_PER_DET=20


import numpy
import pickle

def lookUpDet(allDets, det):
    # change this if we start looking up dets by Dia Id instead of line number!
    return allDets[det]


def readAllDets(mitiFile):
    #this will also have to change if we start looking up dets by Dia Id instead of line number!
    allDets = []
    line = mitiFile.readline()

    while line.strip() != "":
        det = Detection()
        det.fromMitiString(line)
        allDets.append(det)
        line = mitiFile.readline()

    return allDets


class Track(object):
    def __init__(self, trackId, dias, isKnownToBeTrue=None):
        self.trackId = trackId
        self.dias = dias
        self._isTrue = isKnownToBeTrue

    def isTrue(self, dets=None):        
        if self._isTrue != None:
            return self._isTrue
        else:
            if dets == None:
                raise Exception("Can't tell if this is a true track unless it can see the dets, \
                                       or you specified at __init__-time.")
            
            myName = lookUpDet(dets, self.dias[0]).objId
            if myName == 'NS':
                self._isTrue = False
                return False

            for dia in self.dias:
                if lookUpDet(dets, dia).objId != myName:
                    # this is inconsistent, we are not a true track
                    self._isTrue = False
                    return False

            # we didn't find any inconsistency, we are a true track
            self._isTrue = True
            return True

    def isSubset(self, otherTrack):
        myDiasSet = set(self.dias)
        otherTrackDiasSet = set(otherTrack.dias)
        return myDiasSet.issubset(otherTrackDiasSet)
            


class Detection(object):
    def __init__(self, ra=None, dec=None, mjd=None, mag=None, objId=None):
        self.ra = ra
        self.dec = dec
        self.mjd = mjd
        self.mag = mag
        self.objId = objId

    def fromMitiString(self, mitiString):
        items = mitiString.strip().split()
        self.mjd, self.ra, self.dec, self.mag = map(float, items[1:5])
        self.objId = items[6]





class TrackAndRms(object):
    def __init__(self, track, rms):
        self.track = track
        self.rms = rms

        


class TrackQualityBuffer(object):
    """ holds n "best" tracks where n is settable"""

    def __init__(self, nItems):
        self.nItems = nItems
        self.items = [] 

    def add(self, track, rms):
        done = False
        # check if this new track is a superset of any track we hold,
        # and if so, replace the other one.
        for i in range(len(self.items)):
            if self.items[i].track.isSubset(track):
                self.items[i] = TrackAndRms(track=track, rms=rms)
                done = True
        if not done:
            if len(self.items) < self.nItems:
                self.items.append(TrackAndRms(track=track, rms=rms))
            else:
                # find highest-RMS track
                highestRms = None
                highestRmser = None
                for i in range(len(self.items)):
                    if highestRms == None:
                        highestRms=self.items[i].rms
                        highestRmser=i
                    else:
                        rms = self.items[i].rms
                        if rms > highestRms:
                            highestRmser=i
                            highestRms=rms
                #consider replacing the highest-RMS track
                if highestRms > rms:
                    self.items[highestRmser] = TrackAndRms(track=track, rms=rms)

    def getContents(self):
        return self.items







def angularDistance(a, b):
    """ return distance between a and b, where a and b are angles in degrees. """
    while abs(a - b) > 180:
        if a > b:
            b += 360.
        else:
            a += 360.
    return a - b

def convertToStandardDegrees(angle):
    while angle > 360.:
        angle -= 360.
    while angle < 0.:
        angle += 360.
    return angle

def greatCircleDistance(RA0, Dec0, RA1, Dec1):
    """
    return the great-circle distance between two points on the sky,
    uses haversine formula
    """
    deg_to_rad = numpy.pi / 180.
    rad_to_deg = 180. / numpy.pi

    RADist = angularDistance(RA0, RA1);
    DecDist = angularDistance(Dec0, Dec1);    
    #convert all factors to radians
    RADist = deg_to_rad*convertToStandardDegrees(RADist);
    DecDist = deg_to_rad*convertToStandardDegrees(DecDist);
    Dec0 = deg_to_rad*convertToStandardDegrees(Dec0);
    Dec1 = deg_to_rad*convertToStandardDegrees(Dec1);
    r = 2*numpy.arcsin(numpy.sqrt( (numpy.sin(DecDist/2.))**2 + 
                                 numpy.cos(Dec0)*numpy.cos(Dec1)*(numpy.sin(RADist/2))**2) );
    #back to degrees
    return rad_to_deg*r;



def makeContiguous(angles):
    """ given a set of angles (say, RAs or Decs of observation) which
    span a fairly short arc but may actually cross the 0/360 line,
    make these contiguous by using negative angles or whatever is
    necessary.  if this set of angles does NOT span a short arc (>180
    deg) expect all hell to break loose."""
    a0 = angles[0]
    output = [a0]
    for angle in angles[1:]:
        while abs(angle - a0) > 180:
            if angle > a0:
                angle -= 360.
            else:
                angle += 360.
        output.append(angle)
    return output
    
def getRmsForTrack(dets, lineNum):
    t0 = min(map(lambda x: x.mjd, dets))
    ras = []
    decls = []
    mjds = []
    for det in dets:
        ras.append(det.ra)
        decls.append(det.dec)
        mjds.append(det.mjd - t0)
    ras = makeContiguous(ras)
    decls = makeContiguous(decls)
    ras = numpy.array(ras)
    decls = numpy.array(decls)
    mjds = numpy.array(mjds)

    raFunc, raRes, rank, svd, rcond = numpy.polyfit(mjds, ras, 2, full=True)
    decFunc, decRes, rank, svd, rcond = numpy.polyfit(mjds, decls, 2, full=True)
    raFunc = numpy.poly1d(raFunc)
    decFunc = numpy.poly1d(decFunc)

    #now get the euclidean distance between predicted and observed for each point
    netSqDist = 0.0
    dists = []
    for i in range(len(mjds)):
        predRa = raFunc(mjds[i])
        predDec = decFunc(mjds[i])
        dist = greatCircleDistance(predRa, predDec, ras[i], decls[i])
        dists.append(dist)
        if (dist > .1):
            print "Unexpected weirdness at line number %i, diasource had angular distance of %f from best-fit curve prediction" % (lineNum, dist)
            print "Predicted RA, Dec were ", predRa, predDec
            print "observed RA, Dec were ", ras[i], decls[i]
            print "all RAs were ", ras
            print "all decs were ", decls
        sqDist = dist**2
        #print "got euclidean distance was ", sqDist
        netSqDist += sqDist

    rms = numpy.sqrt(netSqDist / len(dets))
    if (rms > .1):
        print "Unexpected weirdness at line number %i, RMS error was %f " % (lineNum, rms)
    return rms, raRes, decRes, dists
 





if __name__ == "__main__":
    import sys
    import os.path
    trackletsMitiFile = file(sys.argv[1],'r')
    tracksByIndicesFile = file(sys.argv[2],'r')
    if (os.path.exists(sys.argv[3])):
        raise Exception("Don't want to open file %s for writing, it appears to exist already")
    tracksToKeepOutFile = file(sys.argv[3], 'w')

    if (os.path.exists(sys.argv[4])):
        raise Exception("Don't want to open file %s for writing, it appears to exist already")
    trackRmsOutFile = file(sys.argv[4], 'w')
    trackRmsOutFile.write("#lineNum isTrueTrack trackLength trackRms dec0\n")

    allDets = readAllDets(trackletsMitiFile)

    findableObjectsInputData = set()
    
    detToBestTracks = {}
    trackLine = tracksByIndicesFile.readline()
    lineNum = 0
    #for each track, see if we keep it
    while trackLine != "":

        #get track detections to calculate RMS
        trackDetsIndices = map(int, trackLine.split())
        trackDets = [ lookUpDet(allDets, x) for x in trackDetsIndices ]
        trackRms, trackRaRms, trackDecRms, perDetDists = getRmsForTrack(trackDets, lineNum)
        thisTrack = Track(lineNum, trackDetsIndices)
        # if this was a true track, add this object to the findable objects set
        isTrue = False
        if thisTrack.isTrue(allDets):
            isTrue = True
            objId = trackDets[0].objId
            findableObjectsInputData.add(objId)

            
        trackRmsOutFile.write("%d %r %d %1.12f %1.12f PER_DET_ERRORS:" % (lineNum, isTrue, len(thisTrack.dias), trackRms, lookUpDet(allDets, thisTrack.dias[0]).dec))
        for dist in perDetDists:
            trackRmsOutFile.write(" %1.12f" % dist)
        trackRmsOutFile.write("\n")

        for det in trackDets:
            if not detToBestTracks.has_key(det):
                detToBestTracks[det] = TrackQualityBuffer(TRACKS_PER_DET)
            detToBestTracks[det].add(track=thisTrack, rms=trackRms)

        trackLine = tracksByIndicesFile.readline()
        if lineNum % 100 == 0:
            print "Just processed track ", lineNum

        lineNum += 1

    print "Done searching through tracks! Read %d tracks." % (lineNum)
    print "There were initially %d findable objects." % (len(findableObjectsInputData))
    
    pickle.dump(detToBestTracks, tracksToKeepOutFile)

    print "Finished dumping dict of det -> TrackQualityBuffer."

    # do some analysis

    tracksSaved = set()
    findableObjectsAfterFilter = set()
    for dets in detToBestTracks.keys():
        tracksAndRmses = detToBestTracks[dets].getContents()        
        for myTrackAndRms in tracksAndRmses:            
            myTrack = myTrackAndRms.track
            tracksSaved.add(myTrack.trackId)
            if myTrack.isTrue():
                objId = lookUpDet(allDets, myTrack.dias[0]).objId
                findableObjectsAfterFilter.add(objId)

    print "Reduced data to %d tracks" % len(tracksSaved)
    print "there were %d findable objects in those tracks" % len(findableObjectsAfterFilter)

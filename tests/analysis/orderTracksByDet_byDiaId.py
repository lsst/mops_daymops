#!/usr/bin/env python

""" jmyers aug 30 2010


Takes the infile and outfile from C++ linkTracklets or a converted
output for C linkTracklets (i.e.: tracks as sets of **dia Ids***, not
line numbers) and then finds the tracks associated with each
detection, ordered by RMS value.  Output some stats about how many
tracks were in the top N for various N.

This version ALSO outputs data about the "approximated" quadratics
(calculated as [last tracklet vel - first tracklet vel] / dt) and the
RMS of detections from THAT approximation.

"""


import numpy
import pickle

def lookUpDet(allDets, det):
    # change this if we start looking up dets by (Dia Id instead of line number!) vice versa really    
    return allDets[det]


def readAllDets(mitiFile):
    #this will also have to change!
    allDets = {}
    line = mitiFile.readline()

    while line.strip() != "":
        det = Detection()
        det.fromMitiString(line)
        allDets[det.diaId] = det
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
    def __init__(self, diaId=None, ra=None, dec=None, mjd=None, mag=None, objId=None):
        self.ra = ra
        self.dec = dec
        self.mjd = mjd
        self.mag = mag
        self.objId = objId

    def fromMitiString(self, mitiString):
        items = mitiString.strip().split()
        self.diaId = int(items[0])
        self.mjd, self.ra, self.dec, self.mag = map(float, items[1:5])
        self.objId = items[6]





class TrackAndRms(object):
    def __init__(self, track, rms):
        self.track = track
        self.rms = rms

        


class TrackQualityBuffer(object):
    """ holds non-subset tracks, ordered by RMS """

    def __init__(self):
        self.items = {}
        self.firstTrueTrackIndex = None

    def add(self, track, rms):
        done = False
        # check if this new track is a superset of any track we hold,
        # and if so, replace the other one.
        if self.items.has_key(rms):
            print "WARNING: redundant RMS of %1.12f" % rms
            print "tracks are: "
            print self.items[rms].dias
            print track.dias
            raise Exception("Multiple tracks have same RMS: %1.12f") % rms
        for key in self.items.keys():
            if self.items[key].track.isSubset(track):
                self.items[key] = TrackAndRms(track=track, rms=rms)
                done = True
        if not done:
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

    def firstTrueTrackRank(self):
        """ return the number of items which have RMS lower than the first true track"""
        count = 0
        for key in sorted(self.items.keys()):
            if items[key].isTrue():
                return count
            count += 1
        #no true tracks at all!
        return None
        

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



def getNetSqDist(ras, decls, mjds, raFunc, decFunc, isApproximate=False):
    netSqDist = 0.0
    dists = []
    for i in range(len(mjds)):
        predRa = raFunc(mjds[i])
        predDec = decFunc(mjds[i])
        dist = greatCircleDistance(predRa, predDec, ras[i], decls[i])
        dists.append(dist)
        if not isApproximate and (dist > .1):
            print "Unexpected weirdness at line number %i, diasource had angular distance of %f from best-fit curve prediction" % (lineNum, dist)
            print "Predicted RA, Dec were ", predRa, predDec
            print "observed RA, Dec were ", ras[i], decls[i]
            print "all RAs were ", ras
            print "all decs were ", decls
        sqDist = dist**2
        #print "got euclidean distance was ", sqDist
        netSqDist += sqDist
    return netSqDist, dists


def getRmsForTrack(dets, lineNum):
    t0 = min(map(lambda x: x.mjd, dets))
    ras = []
    decls = []
    mjds = []
    for det in sorted(dets, key=lambda d: d.mjd):
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
    preciseRaAccel = raFunc[0]
    preciseDecAccel = decFunc[0]

    #now get the euclidean distance between predicted and observed for each point
    netSqDist, dists = getNetSqDist(ras, decls, mjds, raFunc, decFunc)

    rms = numpy.sqrt(netSqDist / len(dets))

    if (rms > .1):
        print "Unexpected weirdness at line number %i, RMS error from best-fit curve was %f " % (lineNum, approximatedRms)

    # also in this version: compute APPROXIMATED RMS by taking the
    # APPROXIMATE acceleration (calculated by comparing endpoint
    # tracklet velocities) and 
    firstTrackletRaVel = (ras[1] - ras[0]) / (mjds[1] - mjds[0])
    firstTrackletDecVel = (decls[1] - decls[0]) / (mjds[1] - mjds[0])

    lastTrackletRaVel = (ras[-1] - ras[-2]) / (mjds[-1] - mjds[-2])
    lastTrackletDecVel = (decls[-1] - decls[-2]) / (mjds[-1] - mjds[-2])

    approximatedRaAccel = (lastTrackletRaVel - firstTrackletRaVel) / (mjds[-2] - mjds[0])
    approximatedDecAccel = (lastTrackletDecVel - firstTrackletDecVel) / (mjds[-2] - mjds[0])

    approximatedRaFunc = numpy.poly1d([ approximatedRaAccel, firstTrackletRaVel, ras[0]])
    approximatedDecFunc = numpy.poly1d([approximatedDecAccel, firstTrackletDecVel, decls[0]])

    netSqDist, approximatedDists = getNetSqDist(ras, decls, mjds, approximatedRaFunc, approximatedDecFunc, isApproximate=True)
    approximatedRms = numpy.sqrt(netSqDist / len(dets))

    return rms, preciseRaAccel, preciseDecAccel, raRes, decRes, approximatedRms, approximatedRaAccel, approximatedDecAccel, dists, approximatedDists
 





if __name__ == "__main__":
    import sys
    import os.path
    detectionsMitiFile = file(sys.argv[1],'r')
    tracksByDiaIdFile = file(sys.argv[2],'r')
    if (os.path.exists(sys.argv[3])):
        raise Exception("Don't want to open file %s for writing, it appears to exist already")
    tracksToKeepOutFile = file(sys.argv[3], 'w')

    if (os.path.exists(sys.argv[4])):
        raise Exception("Don't want to open file %s for writing, it appears to exist already")
    trackRmsOutFile = file(sys.argv[4], 'w')
    trackRmsOutFile.write("#lineNum isTrueTrack trackLength trackRms dec0 raAccel decAccel approxTrackRms approxRaAccel approxDecAccel PER_DET_ERRORS APPROX_PER_DET_ERRORS\n")

    allDets = readAllDets(detectionsMitiFile)

    findableObjectsInputData = set()
    
    detToBestTracks = {}
    trackLine = tracksByDiaIdFile.readline()
    lineNum = 0
    #for each track, see if we keep it
    while trackLine != "":

        #get track detections to calculate RMS
        trackDetsIndices = map(int, trackLine.split())
        trackDets = [ lookUpDet(allDets, x) for x in trackDetsIndices ]

        #old: trackRms, trackRaRms, trackDecRms, perDetDists = getRmsForTrack(trackDets, lineNum)

        [trackRms, preciseRaAccel, preciseDecAccel, trackRaRms, trackDecRms, 
         approximatedRms, approximatedRaAccel, approximatedDecAccel, perDetDists, approximatedDists] =  getRmsForTrack(trackDets, lineNum)

        thisTrack = Track(lineNum, trackDetsIndices)
        # if this was a true track, add this object to the findable objects set
        isTrue = False
        if thisTrack.isTrue(allDets):
            isTrue = True
            objId = trackDets[0].objId
            findableObjectsInputData.add(objId)

            
        trackRmsOutFile.write("%d %r %d %1.12f %1.12f %1.12f %1.12f %1.12f %1.12f %1.12f PER_DET_ERRORS:" % \
                                  (lineNum, isTrue, len(thisTrack.dias), trackRms, lookUpDet(allDets, thisTrack.dias[0]).dec, \
                                       preciseRaAccel, preciseDecAccel, approximatedRms, approximatedRaAccel, approximatedDecAccel))
        for dist in perDetDists:
            trackRmsOutFile.write(" %1.12f" % dist)

        trackRmsOutFile.write(" APPROX_PER_DET_ERRORS: ") 

        for dist in approximatedDists:
            trackRmsOutFile.write(" %1.12f" % dist)

        trackRmsOutFile.write("\n")


        for det in trackDets:
            if not detToBestTracks.has_key(det):
                detToBestTracks[det] = TrackQualityBuffer()
            detToBestTracks[det].add(track=thisTrack, rms=trackRms)

        trackLine = tracksByDiaIdFile.readline()
        if lineNum % 100 == 0:
            print "Just processed track ", lineNum

        lineNum += 1

    print "Done searching through tracks! Read %d tracks." % (lineNum)
    print "There were initially %d findable objects." % (len(findableObjectsInputData))
    
    pickle.dump(detToBestTracks, tracksToKeepOutFile)

    print "Finished dumping dict of det -> TrackQualityBuffer."
    print "doing some analysis."

    # do some analysis
    trueTrackRankHistogram = {}

    tracksSaved = set()    
    findableObjectsAfterFilter = set()

    for dets in detToBestTracks.keys():
        firstTrueTrack = detToBestTracks[dets].getFirstTrueTrackRank()
        if firstTrueTrack != None:
            if trueTrackRankHistogram.has_key(firstTrueTrack):
                trueTrackRankHistogram[firstTrueTrack] += 1
            else:
                trueTrackRankHistogram[firstTrueTrack] = 1

        tracksAndRmses = detToBestTracks[dets].getContents()        
        for myTrackAndRms in tracksAndRmses:            
            myTrack = myTrackAndRms.track
            tracksSaved.add(myTrack.trackId)
            if myTrack.isTrue():
                objId = lookUpDet(allDets, myTrack.dias[0]).objId
                findableObjectsAfterFilter.add(objId)

    print "Reduced data to %d tracks via subset removal" % len(tracksSaved)
    print "there were %d findable objects in those tracks" % len(findableObjectsAfterFilter)

    for rank in sorted(firstTrueTrackRankHistogram.keys()):
        print "%d\t dets had the true track as their #%d best track" % (firstTrueTrackRankHistogram[rank], rank + 1)

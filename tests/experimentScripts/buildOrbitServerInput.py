#!/usr/bin/env python

"""

jmyers

Build input for orbit_server.x.  Takes as input a set of detections in
MITI format and a set of tracks expressed as DiaSourceIds.  outputs a
variety of files with a given prefix, as orbit_server.x enjoys.

Unlike the PS version, we do not require that a track be made of
tracklets but merely of a set of detections.  (Kubica's linkTracklets
will allow this to happen, as will ours, though it is unlikely except
in fairly dense observations.)  To allow this, we actually write
"singleton tracklets" (i.e., single detections) to the ".in.tracklets"
file handed to orbit_server.  This actually appears to work just fine!

output files

prefix.in.manifest : lists the names of the three files we are writing. Not sure why we need this!
prefix.in.tracklets : holds the tracklets (really, the detections) comprising the tracks.
prefix.in.request : holds the tracks, and the flag requesting IOD be performed on them

"""

ASSUMED_OBS_TYPE='O'
ASSUMED_FILTER='r'
ASSUMED_OBSERVATORY=866
ASSUMED_RMS_RA=0.3
ASSUMED_RMS_DEC=0.3
ASSUMED_RMS_MAG=0.1
ASSUMED_S2N=18.0

DEBUG=True

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




def writeTrackletsFile(detsFile, outFile):
    """ translate from MITI detections to orbit_server.x input format
    'tracklets.' we don't really write tracklets at all -
    orbit_server.x will actually allow 'singleton' tracklets of size
    one, which are basically just detections.  The 'tracklets' we
    write have the same IDs as the DiaIds on the detections from the
    input file. easy peasy!
    
    TBD: right now we use constants for OBS_TYPE, FILTER, OBSERVATORY,
    RMS_RA, RMS_DEC, RMS_MAG and S2N. In the future, find these values
    for real (especially if we really need them)
    """
    
    outFile.write("!!OID TIME OBS_TYPE RA DEC APPMAG FILTER OBSERVATORY RMS_RA RMS_DEC RMS_MAG S2N Secret_name\n")

    detsLine = detsFile.readline()
    while detsLine != "":
        det = Detection()
        det.fromMitiString(detsLine)
        outFile.write("%d %5.10f %s %3.12f %3.12f %3.12f %s %s %3.12f %3.12f %3.12f %3.12f %s\n" \
                          % \
                      (det.diaId, det.mjd, ASSUMED_OBS_TYPE, det.ra, det.dec, det.mag, ASSUMED_FILTER,
                       ASSUMED_OBSERVATORY, ASSUMED_RMS_RA, ASSUMED_RMS_DEC, ASSUMED_RMS_MAG,
                       ASSUMED_S2N, det.objId))
        detsLine = detsFile.readline()


def writeRequest(inTracks, requestFile):
    """ write the .in.request file given a set of DIA IDs which comprise the tracks."""
    requestFile.write("!!ID_OID NID TRACKLET_OIDs OP_CODE N_OBS N_SOLUTIONS N_NIGHTS ARC_TYPE NO_RADAR PARAM(4)\n")
    trackLine = inTracks.readline()

    count = 0
    while trackLine != "":        
        diaIds = map(int, trackLine.split())
        trackName = '='.join(map(str, diaIds)) # this is what PS did so just run with it i guess
        requestFile.write("%s" % trackName)
        requestFile.write(" %d" % len(diaIds))
        for dia in diaIds:
            requestFile.write(" %d" % dia)
        requestFile.write(" REQUEST_PRELIM %d 0 0 0 0 0.0 0.0 0.0 0.0\n" % len(diaIds))
        trackLine = inTracks.readline()
        count += 1
        if DEBUG and count > 1000:
            return



if __name__=="__main__":
    import sys
    inDets = file(sys.argv[1],'r')
    inTracks= file(sys.argv[2],'r')
    outPrefix=sys.argv[3]

    
    manifestName = outPrefix + '.in.manifest'
    trackletsName = outPrefix + '.in.tracklet'
    requestName = outPrefix + '.in.request'
    manifest = file(manifestName, 'w')

    # write the manifest file. too connected and too simple to really put in its own function.
    for name in [manifestName, trackletsName, requestName]:
        manifest.write(name + '\n')

    trackletsOut = file(trackletsName,'w')
    writeTrackletsFile(inDets, trackletsOut) 

    requestFile = file(requestName,'w')
    writeRequest(inTracks, requestFile)
    

#!/usr/bin/env python

"""

jmyers


Build input(s) for orbit_server.x.  This script takes as input a set
of tracks expressed as DiaSourceIds; the DiaSources are looked up from
a database.  It then outputs a variety of files with a given prefix,
as orbit_server.x enjoys.  The outputs are sets of files to be used
for orbit_server.x.  Because orbit_server.x has a hard-coded maximum
number of items in its input files, it is usually necessary to write
many sets of files.


output files

prefix_num.in.manifest : lists the names of the three files we are writing. Not sure why we need this!
prefix_num.in.tracklets : holds the tracklets (really, the detections) comprising the tracks.
prefix_num.in.request : holds the tracks, and the flag requesting IOD be performed on them

orbit_server.x takes a maximum of 500000 items in the .in.tracklets file.

It also appears to take a maximum of about 24000 items in the .in.request file.

We will create as many sets of files as needed.

Unlike the PS version, we do not require that a track be made of
tracklets but merely of a set of detections.  (Kubica's linkTracklets
will allow this to happen, as will ours, though it is unlikely except
in fairly dense observations.)  To allow this, we actually write
"singleton tracklets" (i.e., single detections) to the ".in.tracklets"
file handed to orbit_server.  This actually appears to work just fine
- I tested by taking one of the PS examples and renaming the tracklets
and changing the contents of the tracks and got the same results out.


"""

ASSUMED_OBS_TYPE='O'
ASSUMED_FILTER='r'
ASSUMED_OBSERVATORY=807
ASSUMED_RMS_MAG=1.
OPSIM_DB='opsim_3_61'
OPSIM_TABLE='output_opsim3_61'
#DIAS_DB='mops_noDeepAstromError'
#DIAS_TABLE='fullerDiaSource'
DIAS_DB='fullSkyOneMonth'
DIAS_TABLE='mopsDetections'

MAX_TRACKLETS_PER_FILE=50000
MAX_TRACKS_PER_FILE=24000

# set MAX_TRACKS to -1 to allow infinite tracks in output. otherwise
# set to a fixed number to do a smaller, scaled run.  handy for
# debugging.
MAX_TRACKS=-1

import add_astrometric_noise as astrom
import useful_input as ui
from lsst.daf.base.baseLib import DateTime



def writeOptFile(optFileName):
    s = """! option file for OrbFit Server Suite   
orbsrv.  

! oss1run1 
! slow iteration of orbit_server with input tracks from MOPS code 
     .force_difcor=.T.       ! compute least sq. orbit  
!     .inpdir='.'            ! input directory (def .)
!     .outdir='.'            ! output directory (def .)
     .cooy='COM'             ! output cometary orbital elements (COM)
! logical controls of prelim_orbit  
   .smart_ecclim_cont=0.1d0  ! eta max for hyperbolics, deg/d (def. 0)
   .smart_gaussmap=.T.       ! use elongation to select Gauss map (def F)
! control parameters for prelim_orbit  
     .prelim_rms=  1000.d0 ! output prelim if RMS<this (def 100)  
!

! control parameters for fit_control
     .rms_control=1.8d0     ! control on astrometric RMS (def 1.8)
     .rmsh_control=0.4d0    ! control on photometric RMS (def 1.5)
     .bias_control=4.d0     ! control on bias of the fit (def 2.5)
     .span_control=4.d0     ! control on linear trends (def 2.5)
     .curv_control=4.d0     ! control on residual curvature (def 2.5)
     .zsign_control=4.d0    ! control on residual third deriv. (def 2.5)
"""
    out = file(optFileName, 'w')
    out.write(s)



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





def taiToUtc(taiTime):
    d = DateTime(taiTime, DateTime.TAI)
    return d.mjd(DateTime.UTC)


def lookUpDet(dbCursor, diaId):
    """ fetches the following from DB, given a DiaID, doing a little calculation as
    needed: utcMjd, ra, dec, mag, astromErr, snr, groundTruthId
    """
        
    sql = """ SELECT taiMidPoint, ra, decl, mag, opSimId, snr, ssmId
                  FROM %s.%s WHERE diaSourceId=%d ; """ % \
        (DIAS_DB, DIAS_TABLE, diaId)
    res = ui.sqlQuery(dbCursor, sql, verbose=False)
    #print res
    [taiTime, ra, dec, mag, opSimId, snr, groundTruthId] = res[0]
    #print taiTime, ra, dec, mag, opSimId, snr, groundTruthId
    sql = "select 5sigma_ps, seeing from %s.%s where obsHistID=%d" %(OPSIM_DB, OPSIM_TABLE, opSimId)
    [sigma5_ps, seeing] = ui.sqlQuery(cursor, sql, verbose=False)[0]
    
    astromErr = astrom.calcAstrometricError(mag, sigma5_ps, seeing)

    #print "callint taiToUtc(", taiTime, ")"
    #utcMjd = taiToUtc(taiTime)
    #print " calling DateTime(", taiTime, ", DateTime.TAI)"
    d = DateTime(taiTime, DateTime.MJD, DateTime.TAI)
    utcMjd = d.mjd(DateTime.UTC)
    #print "got tai ", taiTime, " = utc ", utcMjd
    return utcMjd, ra, dec, mag, astromErr, snr, groundTruthId




def writeDiasToDesTrackletsFile(diasToWrite, outFile, dbCursor):
    """ reads a set of dias from the DB and writes them to the
    "tracklets" file. we don't really write tracklets at all -
    orbit_server.x will actually allow 'singleton' tracklets of size
    one, which are basically just detections.  The 'tracklets' we
    write have the same IDs as the DiaIds on the detections from the
    input file. easy peasy!
    """
    for diaId in diasToWrite:
        utcMjd, ra, dec, mag, astromErr, snr, groundTruthId = lookUpDet(dbCursor, diaId)

        #TBD: RA astrometric error is not == astromErr, it's something else - cos(radians(dec))* astromErr?
        #   add this in later. For now on the ecliptic it won't really matter.

        outFile.write("%d %5.10f %s %3.12f %3.12f %3.12f %s %s %3.12f %3.12f %3.12f %3.12f %s\n" \
                          % \
                      (diaId, utcMjd, ASSUMED_OBS_TYPE, ra, dec, mag, ASSUMED_FILTER,
                       ASSUMED_OBSERVATORY, astromErr, astromErr, ASSUMED_RMS_MAG,
                       snr, groundTruthId))


def writeTrackToRequestFile(diaIds, requestFile, trackName):
    """ take a track ( as a set of diaIds ) and write it to request file."""

    t = Track(trackName, diaIds)

    #if this track has too many dias, then we need to subselect. max
    #is 18 according to orbit_server.x

    if len(diaIds) > 18:
        half = len(diaIds)/2
        newDias = diaIds[:8] + diaIds[half:half+2] + diaIds[-8:]
        diaIds = newDias
        if len(diaIds) > 18: 
            print "whoops! programmer error in writeTracksToRequestFile"
    #outputTrack = True
    #if TRUE_TRACKS_ONLY == True:
    #    outputTrack = t.isTrue(allDets)
    #if FALSE_TRACKS_ONLY== True:
    #    outputTrack = not t.isTrue(allDets)
    #if outputTrack:
    requestFile.write("%s" % trackName)
    requestFile.write(" %d" % len(diaIds))
    for dia in diaIds:
        requestFile.write(" %d" % dia)
    requestFile.write(" REQUEST_PRELIM %d 0 0 0 0 0.0 0.0 0.0 0.0\n" % len(diaIds))




def writeManifestFile(manifestName, trackletsName, requestName):
    manifest = file(manifestName,'w')
    for name in [manifestName, trackletsName, requestName]:
        manifest.write(name + '\n')




def createNewFiles(outPrefix, outputNumber):
    outNumS = str(outputNumber)
    requestFileName = outPrefix + "_" + outNumS + ".in.request"
    trackletsFileName = outPrefix + "_" + outNumS + ".in.tracklet"
    manifestName = outPrefix + "_" + outNumS + ".in.manifest"
    writeManifestFile(manifestName, trackletsFileName, requestFileName)

    optFileName = outPrefix + "_" + outNumS + ".opt"
    writeOptFile(optFileName)

    requestFile = file(requestFileName,'w')
    trackletsFile = file(trackletsFileName,'w')
    trackletsFile.write("!!OID TIME OBS_TYPE RA DEC APPMAG FILTER OBSERVATORY RMS_RA RMS_DEC RMS_MAG S2N Secret_name\n")
    requestFile.write("!!ID_OID NID TRACKLET_OIDs OP_CODE N_OBS N_SOLUTIONS N_NIGHTS ARC_TYPE NO_RADAR PARAM(4)\n")

    return requestFile, trackletsFile



def writeOrbitServerInputFiles(inTracksFile, outPrefix, cursor, maxTrackletsPerFile, maxTracksPerFile):
    """writes sets of input files for orbit_server.x. create as many sets of output files as needed"""
    curFileSetNum = 0
    totalTracksWritten = 0
    requestFile, trackletsFile = createNewFiles(outPrefix, curFileSetNum)
    
    nextTrackLine = inTracksFile.readline()

    diasNeededForThisOutfile = set()
    numTracksWrittenToThisOutfile = 0

    while nextTrackLine != "":

        nextTrack = map(int, nextTrackLine.split())
        # we need to know the next track in order to see if it will
        # require us to add too many detections to an outfile.
        readaheadTrackLine = inTracksFile.readline()
        readaheadTrack = map(int, readaheadTrackLine.split())

        # check if we're done with this set of files for any reason.
        # possible reasons:
        #  - we are out of input (done)
        #  - we can't fit any more detections in an outfile.
        #  - we can't fit any more *tracks* in an outfile.
        if readaheadTrackLine == "" or \
                len(diasNeededForThisOutfile) + len(readaheadTrack) > maxTrackletsPerFile or \
                numTracksWrittenToThisOutfile == maxTracksPerFile or \
                (MAX_TRACKS > 0 and totalTracksWritten == MAX_TRACKS):

            # if we're done with this set of files, start a new set of files.
            writeDiasToDesTrackletsFile(diasNeededForThisOutfile, trackletsFile, cursor)
            requestFile.close()
            trackletsFile.close()
            print "Successfully wrote output files for %s" % (outPrefix + "_" + str(curFileSetNum))
            curFileSetNum += 1
            requestFile, trackletsFile = createNewFiles(outPrefix, curFileSetNum)
            diasNeededForThisOutfile = set()
            numTracksWrittenToThisOutfile = 0
            

        # write this track to the current outfile; add its dias to the set needed for the associated
        # "tracklets" file (which is really a set of detections)        
        for dia in nextTrack: 
            diasNeededForThisOutfile.add(dia)
        writeTrackToRequestFile(nextTrack, requestFile, str(totalTracksWritten))
        totalTracksWritten += 1
        numTracksWrittenToThisOutfile += 1

        # collect more data for the current file.
        nextTrackLine = readaheadTrackLine

    requestFile.close()
    trackletsFile.close()


if __name__=="__main__":
    import sys
    inTracks= file(sys.argv[1],'r')
    outPrefix=sys.argv[2]

    conn,cursor = ui.sqlConnect(hostname='localhost', username='jmyers', passwdname='jmyers', dbname='opsim_3_61')
    
    writeOrbitServerInputFiles(inTracks, outPrefix, cursor, MAX_TRACKLETS_PER_FILE, MAX_TRACKS_PER_FILE)
    print "done writing output files."

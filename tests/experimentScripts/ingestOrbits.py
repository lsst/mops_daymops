#!/usr/bin/env python

""" jmyers march 19, 2012

this is a script for reading the output of orbit_server.x into a DB
table for later analysis/attribution/precovery. 

Say you had a set of tracks e.g. my.tracks built from diaSources dias.miti.  You then used
buildOrbitServerInput.py to create lots of files for orbit_server.
Then you used runOrbitServer.sh and got lots of .out.orbit files.

Give this script the following:

my.tracks dias.miti path/to/\*.orbit.out 

Remember to escape the \* (or whatever regexp you use) - we want
Python to expand it, not the shell.

"""



import mopsDatabases
import glob, sys



def writeOutput(trackId, trackDias, orbit, dias):
    """ this is called once for every track/orbit combination, and
    also once with orbit=None for every failed IOD

    trackId is really only used for debug purposes.
    """
    # each d in dias is [diaId, mjd, ssmId]
    tds = map(lambda x: dias[x], trackDias)
    mjds = map(lambda x: x[1], tds)
    #print tds
    ssmIds = map(lambda x: x[2], tds)
    if len(set(ssmIds)) > 1:
        ssmId = -1
    else:
        ssmId = ssmIds[0]
    dt = max(mjds) - min(mjds)
    t0 = min(mjds)
    if orbit==None:
        print "track %d: objId = %d, dt=%f, t0=%f orbit: %s" % (trackId,ssmId, dt, t0, "IOD_FAILED")

    else:
        print "track %d: objId = %d, dt=%f, t0=%f orbit: %s" % (trackId,ssmId, dt, t0, orbit.rstrip())
    




def cmpFilenames(n1, n2):
    """ standard string comparison thinks that my_20.out.orbit >
    my_100.out.orbit.  This is a different comparator which knows that
    longer filenames correspond to larger numbers.  Obviously this
    only works for filenames which are identical except for the
    numbered part."""
    if len(n1) > len(n2):
        return 1
    elif len(n1) < len(n2):
        return -1
    else:
        # strings of the same length can be compared normally
        return cmp(n1, n2)



def getTrackStats(trackDias, curs):
    """ return underlying object ID (or -1 if a false track) and the
    start time, dt of the track, using DB"""
    

def getNextOrbit(orbitsFiles):
    """helper function for ingestOrbits. Find the next non-empty,
    non-header line in the set of files, which MUST be sorted. If
    there is no such line return empty string."""
    f = orbitsFiles[0]
    l = f.readline()
    while True:
        if l != "":
            if l[0] != "!":
                # first case: we try to read a line from the current
                # file, and there is one, and it is not a header, so
                # return it.
                orbitsFiles[0] = f
                return l, orbitsFiles
            else:
                # this is a header. try again.
                #print "Saw a header line"
                l = f.readline()
        else:
            # we just hit the end of a file.  throw away this file,
            # then get the next one
            orbitsFiles = orbitsFiles[1:]
            if orbitsFiles == []:
                print "DONE. No more orbits files left."
                sys.exit(0)
            else:
                #print "moving on to next file, ", orbitsFiles[0], "... there are now ", \
                #    len(orbitsFiles), " files to look in"
                f = orbitsFiles[0]
                l = f.readline()

def ingestOrbits(tracksFile, orbitsFiles, dias):
    """ REQUIRE and ASSUME that orbitsFiles are sorted; that is,
    orbits from tracks 1 thru n are in the first file and n+1 thru n+m
    are in the second file and so on.  We must have at least as many
    input tracks as all the entries in the orbitsFiles."""
    curTrack = 0
    curOrbit, orbitsFiles = getNextOrbit(orbitsFiles)
    curOrbitTrack = int(curOrbit.split()[0])
    trackLine = tracksFile.readline()
    while trackLine != "":
        trackDias = map(int, trackLine.split())
        if curOrbitTrack == curTrack:
            # we got the track corresponding to this orbit.
            #print "   CurTrack = ", curTrack
            #print "   ", curOrbit
            curOrbit, orbitsFiles = getNextOrbit(orbitsFiles)            
            curOrbitTrack = int(curOrbit.split()[0])

            writeOutput(curTrack, trackDias=trackDias, \
                            orbit=curOrbit, dias=dias)
            while curTrack == curOrbitTrack: 
                # sometimes we get multiple orbit solutions for a
                # single track
                writeOutput(curTrack, trackDias=trackDias, \
                                orbit=curOrbit, dias=dias)

                curOrbit, orbitsFiles = getNextOrbit(orbitsFiles)            
                curOrbitTrack = int(curOrbit.split()[0])
        else:
            # IOD failed for this track
            writeOutput(curTrack, trackDias=trackDias, orbit=None, dias=dias)
        trackLine = tracksFile.readline()
        curTrack += 1
    print "DONE.  Saw ", curTrack+1, " tracks"



def readDias(diasDataFile):
    """ reads a dump of dias, which include diaId expMjd ssmId
    obsHistId for every diaSource.  Returns a dict mapping diaId to
    data."""
    print "Reading all diaSources into memory..."
    idToDias = {}
    line = diasDataFile.readline()
    while line != "":
        [diaId, expMjd, ra, decl, mag, obscode, ssmId, e, ee] = \
            line.split()
        diaId = int(diaId)
        expMjd = float(expMjd)
        ssmId = int(ssmId)
        
        idToDias[diaId] = [diaId, expMjd, ssmId]
        
        line = diasDataFile.readline()

    print "Finished reading diaSources."
    return idToDias



if __name__=="__main__":

    tracksFile, diasFile, orbitFiles = sys.argv[1:]
    orbitFiles = glob.glob(orbitFiles)
    dias = readDias(file(diasFile))
    # get the .out.orbit files in order - first prefix_1.out.orbit,
    # then prefix_2.out.orbit, then prefix_100.out.orbit.  normal
    # Python cmp will give these in the wrong order, use a custom
    # comparator.  We do assume that these all share the same prefix    
    # and suffix...
    orbitFiles = sorted(orbitFiles, cmp=cmpFilenames)
    orbitFiles = map(lambda x: file(x,'r'), orbitFiles)
    ingestOrbits(file(tracksFile,'r'), orbitFiles, dias)

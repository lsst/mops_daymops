#!/usr/bin/env python

"""

jmyers oct 15

Bin tracklets into multiple outfiles; write the start_t_range value
for each one.  Create each outfile s.t. it contains data from one
night or at most some fixed number of start images; this will improve
load-balancing per run and hopefully keep us from exhausting our
per-process allocations on abe.

our old script for making input files for linkTracklets
(make15DayWindows.py) was based on the assumption that we would run
linkTracklets sequentially, without a per-process time limit.

On Abe, we have to fit all our work into a known time limit before the
cluster kills our job.  Also, we would like to get better parallelism.


"""
import os.path
import MySQLdb as db
import sys
import glob
import time

import mopsDatabases



# obscode added to MITI files.
FORCED_OBSCODE="807"

# this needs to be larger than the min time between images.
EPSILON=1e-5

MAX_START_IMAGES_PER_RUN=30

TRACKLETS_BY_OBSHIST_SUFFIX=".tracklets.byDiaId"
TRACKLETS_FILE_TO_OBSHIST=lambda x: (int(os.path.basename(x)[:-len(TRACKLETS_BY_OBSHIST_SUFFIX)]))

# if true, will put some metadata about each data set into the start t range files directory.
WRITE_ADDITIONAL_METADATA=True

# jmyers: this seems to work for our data set but test your own first!
MJD_TO_NIGHT_NUM=lambda mjd: int(mjd)


# set to True if dias will fit in memory; they will load at beginning
# of execution and this will allow for MUCH faster lookups later.
# set to False if debugging to avoid a big wait at startup.
PRELOAD_DIAS=True
PRELOAD_DIAS_FROM_FILE=True



WRITE_CPP_INFILES=True
WRITE_C_INFILES=False


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



def obsHistToTrackletsFile(obsHist, obsHistDir):
    return os.path.join(obsHistDir) + str(obsHist) + TRACKLETS_BY_OBSHIST_SUFFIX


def fetchAllDiasFromFile(diasTableOrFile):
    """ return a LIST mapping diaId to diaSource, for all Dias in the
    given file. This fails if DiaIds do not ascend contiguously from 0
    or 1 to their max value"""
    toRet = [0]
    print "Reading DiaSources from file..."
    f = file(diasTableOrFile,'r')
    nLines = 0
    line = f.readline()
    #fill the list with nLines elements
    while line != "":
        nLines += 1
        toRet.append(0)
        line = f.readline()
    f.close()
    print "File contains ", nLines, " elements. Created an ", len(toRet), " list to hold them."
    #now populate the list
    f = file(diasTableOrFile,'r')
    line = f.readline()
    count = 0;
    while line != "":
        diaId, obshistId, ssmId, ra, decl, MJD, mag, snr = line.split()
        diaId, obshistId, ssmId = map(int, [diaId,obshistId, ssmId])
        ra, decl, MJD, mag, snr = map(float, [ra, decl, MJD, mag, snr])
        #d = Detection(diaId=diaId, ra=ra, dec=decl, mjd=MJD, mag=mag, objId=ssmId)
        toRet[diaId] = [diaId, ra, decl, MJD, mag, ssmId]
        line = f.readline()
        count += 1
        if count % 100000 == 0:
            print "Read ", count, " dias so far..."
    print "Done reading diaSources."
    f.close()
    return toRet


def fetchAllDiasFromDb(cursor, diasDb, diasTable):
    """ return a dictionary mapping diaId to diaSource, for all Dias
    in our data set."""
    toRet = {}
    print "Fetching DiaSources from DB into memory, a few at a time..."
    s = """ SELECT diaSourceId, ra, decl, taiMidPoint, mag, ssmId FROM
            %s.%s; """ % (diasDb, diasTable)
    cursor.execute(s)
    print "  got results from DB, converting them to Python objects..."
    for row in cursor:
        [diaId, ra, dec, mjd, mag, objId] = row
        d = Detection(diaId=diaId, ra=ra, dec=dec, mjd=mjd, mag=mag, objId=objId)
        toRet[diaId] = d
        numFetched = len(toRet)
        if numFetched % 100000 == 0:
            print "  Fetched ", numFetched, " dias so far "
            print "  Time is ", time.ctime()
            sys.stdout.flush()

    print "... Done fetching dias."
    return toRet




def lookUpDia(allDias, diaId, cursor=None, diasDb=None, diasTable=None):
    if PRELOAD_DIAS:
        if PRELOAD_DIAS_FROM_FILE:
            [diaId, ra, dec, mjd, mag, objId] = allDias[diaId]
            d = Detection(diaId=diaId, ra=ra, dec=dec, mjd=mjd, mag=mag, objId=objId)
            return d
        else:
            return allDias[diaId]
    else:
        if None in [cursor, diasDb, diasTable]:
            raise Exception("lookUpDia: If dias are not preloaded, you must specify a cursor, dias table and dias DB name")
        s = """ SELECT diaSourceId, ra, decl, taiMidPoint, mag, ssmId FROM 
                %s.%s WHERE diaSourceId=%d; """ % (diasDb, diasTable, diaId)
        cursor.execute(s)
        [diaId, ra, dec, mjd, mag, objId] = cursor.fetchone()
        d = Detection(diaId=diaId, ra=ra, dec=dec, mjd=mjd, mag=mag, objId=objId)
        return d



def getExpMjdAndFieldIdForObsHists(cursor, obsHists):

    """ use the DB to look up all the obsHists in obsHists and make a
    dictionary mapping obsHistId to time of image. Also make a dict
    mapping obsHistId to fieldId. Return both dictionaries in that
    order."""

    s = """ SELECT expMjd, fieldId, obsHistId FROM %s.%s where obsHistId in ( """ \
        % (mopsDatabases.OPSIM_DB, mopsDatabases.OPSIM_TABLE)

    for i in range(len(obsHists)):
        obsHist = obsHists[i]
        s += str(obsHist)
        if i != len(obsHists) - 1:
            s += ", " 
    s += ");"

    cursor.execute(s)
    results = cursor.fetchall()
    
    obsHistToExpMjd = {}
    obsHistToFieldId = {}
    for row in results:
        obsHistToExpMjd[row[2]] = row[0]
        obsHistToFieldId[row[2]] = row[1]
    return obsHistToExpMjd, obsHistToFieldId



def getRemainingObsHistsForNight(nightNumToObsHists, assignedObsHists, nightNum):
    if not nightNumToObsHists.has_key(nightNum):
        return []

    obsHistsTonight = nightNumToObsHists[nightNum]
    remainingObsHistsTonight = filter(lambda x: x not in assignedObsHists, obsHistsTonight)
    return remainingObsHistsTonight
    


def makeLinkTrackletsInfiles(dbcurs, trackletsByObsHistDir, 
                             outputLinkTrackletsInfilesDir, diasDb, diasTableOrFile, 
                             trackingWindowDays):

    # fetch all dias into memory (just use DB). This will be fasted in
    # the long run than doing lots of cross-language/client-server
    # calls to the DB when we need them.
    if PRELOAD_DIAS:
        if PRELOAD_DIAS_FROM_FILE:
            allDias = fetchAllDiasFromFile(diasTableOrFile)
        else:
            allDias = fetchAllDiasFromDb(dbcurs, diasDb, diasTableOrFile)
    else:
        allDias = None

    # fetch all images which generated tracklets into memory (use a
    # glob of the obshist_to_tracklets directory) get MJDs of each
    # image (use DB)
    trackletsFiles = glob.glob(os.path.join(trackletsByObsHistDir, '*' + TRACKLETS_BY_OBSHIST_SUFFIX))
    allObsHists = map(TRACKLETS_FILE_TO_OBSHIST, trackletsFiles)
    obsHistToExpMjd, obsHistToFieldId = getExpMjdAndFieldIdForObsHists(dbcurs, allObsHists)

    # bin images by night number
    allExpMjds = [ obsHistToExpMjd[oh] for oh in obsHistToExpMjd.keys() ]
    allNightNums = map(MJD_TO_NIGHT_NUM, allExpMjds)

    nightNumsToObsHists = {}
    for nn in allNightNums:
        nightNumsToObsHists[nn] = []
    for obsHist in obsHistToExpMjd.keys():
        nightNumsToObsHists[MJD_TO_NIGHT_NUM(obsHistToExpMjd[obsHist])].append(obsHist)

    # sort images per night by exposure time
    for nn in nightNumsToObsHists.keys():
        nightNumsToObsHists[nn] = sorted(nightNumsToObsHists[nn], key=lambda x: obsHistToExpMjd[x])

    curNightNum = min(allNightNums)
    assignedObsHists = []

    # while there are images not in an output data set:
    while (len(assignedObsHists) < len(allObsHists)):        

        # find the next night with images which need to be assigned. 
        remainingObsHistsTonight = getRemainingObsHistsForNight(nightNumsToObsHists, assignedObsHists, curNightNum)
        while remainingObsHistsTonight == []:
            curNightNum += 1
            remainingObsHistsTonight = getRemainingObsHistsForNight(nightNumsToObsHists, assignedObsHists, curNightNum)

        # choose some images and assign them.
        obsHistsThisDataSet = remainingObsHistsTonight[:MAX_START_IMAGES_PER_RUN]
        assignedObsHists += obsHistsThisDataSet
        lastObsHistTime = min( [ obsHistToExpMjd[oh] for oh in obsHistsThisDataSet ] )

        # find all tracklets from the tonight's data set and next N nights, save them to set of support images
        supportObsHists = []
        #first get the ones from tonight
        for oh in remainingObsHistsTonight:
            if oh not in obsHistsThisDataSet:
                supportObsHists.append(oh)
        # now get those from nights n + 1 throught n + trackingWindowDays
        supportNights = range(curNightNum + 1, curNightNum+trackingWindowDays)
        for supportNight in supportNights:
            if nightNumsToObsHists.has_key(supportNight):
                supportObsHists += nightNumsToObsHists[supportNight]

        # if we have ANY support images at all, then go ahead and write output
        if len(supportObsHists) > 0:
            writeOutputFiles(allDias, dbcurs, obsHistsThisDataSet, supportObsHists, obsHistToExpMjd, 
                             obsHistToFieldId, outputLinkTrackletsInfilesDir, trackletsByObsHistDir, 
                             diasDb, diasTableOrFile)




def writeMitiTrackletsToOutFile(mitiOut, mitiCheatSheetFile, allDias, cursor, trackletsFile, 
                                firstTrackletId, diasDb, diasTable):

    """ expects trackletsFile to be an open file full of tracklets in
    the diaId form (spaces between diaIds, newlines between
    tracklets).  mitiOut and mitiCheatSheetFile are expected to be an
    open file as well, where output is written.  curTrackletId is the
    first ID which will be used.  Returns the number of tracklets
    written."""

    tletsWritten = 0

    line = trackletsFile.readline()
    while line != "":
        items = line.split()
        diaIds = map(int, items)
        for diaId in diaIds:
            #ID EPOCH_MJD RA_DEG DEC_DEG MAG OBSCODE OBJECT_NAME LENGTH ANGLE [EXPOSURE_TIME]

            det = lookUpDia(allDias, diaId, cursor=cursor, diasDb=diasDb, diasTable=diasTable)
            mitiOut.write("%d %2.10f %2.10f %2.10f %2.10f %s %s 0.0 0.0\n" % \
                              (firstTrackletId + tletsWritten,
                               det.mjd, det.ra, det.dec, det.mag, FORCED_OBSCODE, 
                               det.objId))
            mitiCheatSheetFile.write("%d\n" % diaId)

        line = trackletsFile.readline()

        tletsWritten += 1

    return tletsWritten

    
def writeDetsIdsFiles(detsOutFile, idsOutFile, allTrackletsFileNames, allDias, cursor,
                      diasDb, diasTable):
    """ detsOut and idsOut are expected to be open files where output is written.

    writes the simpler linkTracklets input as would be used for C++
    linkTracklets; dets are all the diaSources (miti format) for the
    data set and ids are the tracklets (as sets of LINE NUMBERS into
    the dtes file) for all tracklets in the data set."""

    # first figure out what diaIds need to be written.
    allIds = set()
    for trackletsFileName in allTrackletsFileNames:
        trackletsFile = file(trackletsFileName, 'r')
        tletLine = trackletsFile.readline()
        while tletLine != "":
            ids = map(int, tletLine.split())
            for i in ids:
                allIds.add(i)
            tletLine = trackletsFile.readline()
        trackletsFile.close()

    # write dias that are needed.
    lineNum = 0
    diaToLineNum = {}
    for diaId in allIds:            
        det = lookUpDia(allDias, diaId, cursor=cursor, diasDb=diasDb, diasTable=diasTable)
        detsOutFile.write("%d %2.10f %2.10f %2.10f %2.10f %s %s 0.0 0.0\n" % \
                              (diaId,
                               det.mjd, det.ra, det.dec, det.mag, FORCED_OBSCODE, 
                               det.objId))
        diaToLineNum[diaId] = lineNum
        lineNum += 1

    # now translate from diaIds to indices (line nums) into the file.
    for trackletsFileName in allTrackletsFileNames:
        trackletsFile = file(trackletsFileName, 'r')
        tletLine = trackletsFile.readline()
        while tletLine != "":
            ids = map(int, tletLine.split())
            lineNums = map(lambda x: diaToLineNum[x], ids)
            trackletsAsString = " ".join(map(str, lineNums))
            idsOutFile.write("%s\n" % trackletsAsString)
            tletLine = trackletsFile.readline()
        trackletsFile.close()

def countLines(openFile):
    """ count lines in a file. read one line at a time to avoid memory
    issues if the file is large"""
    count = 0
    line = openFile.readline()
    while line != "":
        count += 1
        line = openFile.readline()
    return count



def writeOutputFiles(allDias, cursor, obsHistsThisDataSet, supportObsHists, obsHistToExpMjd, 
                     obsHistToFieldId, outputLinkTrackletsInfilesDir, trackletsByObsHistDir, 
                     diasDb, diasTable):

    """ write tracklets from the following obsHists into outfiles
    with locations specified by constants in this file."""    

    basename = "image_%d_through_%d" % (obsHistsThisDataSet[0], obsHistsThisDataSet[-1])

    print "Writing output file for ", basename, " at ", time.ctime()
    sys.stdout.flush()

    if WRITE_C_INFILES:
        mitiOutName = basename + ".miti"
        mitiCheatSheetName = basename + ".miti.diaIds"
        mitiOut = file(os.path.join(outputLinkTrackletsInfilesDir, mitiOutName),'w')
        mitiCheatSheetFile = file(os.path.join(outputLinkTrackletsInfilesDir, mitiCheatSheetName),'w')

        curTrackletId = 0
        obsHistToNumTlets = {}
        for obsHist in obsHistsThisDataSet + supportObsHists:
            trackletsFileName = obsHistToTrackletsFile(obsHist, trackletsByObsHistDir)
            trackletsFile = file(trackletsFileName, 'r')
            tletsWritten = writeMitiTrackletsToOutFile(mitiOut, mitiCheatSheetFile, allDias, 
                                                       cursor, trackletsFile, curTrackletId,
                                                       diasDb, diasTable)
            curTrackletId += tletsWritten
            trackletsFile.close()
            obsHistToNumTlets[obsHist] = tletsWritten

        mitiOut.close()
        mitiCheatSheetFile.close()

        # turns out C linkTracklets takes start_t_range as a time offset (in days) from the first image. 
        startTRangeOut = os.path.join(outputLinkTrackletsInfilesDir, basename + ".start_t_range")
        startTRangeOutFile = file(startTRangeOut,'w')
        startTRangeOutFile.write("%f"%(obsHistToExpMjd[obsHistsThisDataSet[-1]] - obsHistToExpMjd[obsHistsThisDataSet[0]] + EPSILON))
        startTRangeOutFile.close()
    # end WRITE_C_INFILES
    else:
        # if not writing C files, we still need to count the number of
        # tracklets per image for our stats file.
        
        obsHistToNumTlets = {}
        for obsHist in obsHistsThisDataSet + supportObsHists:
            trackletsFileName = obsHistToTrackletsFile(obsHist, trackletsByObsHistDir)
            trackletsFile = file(trackletsFileName, 'r')
            tlets = countLines(trackletsFile)
            trackletsFile.close()
            obsHistToNumTlets[obsHist] = tlets


    if WRITE_CPP_INFILES:
        # new: write C++ style outputs
        allTrackletsFiles = []
        for obsHist in obsHistsThisDataSet + supportObsHists:
            trackletsFileName = obsHistToTrackletsFile(obsHist, trackletsByObsHistDir)
            allTrackletsFiles.append(trackletsFileName)

        detsOutName = basename + ".dets"
        detsOutFile = file(os.path.join(outputLinkTrackletsInfilesDir, detsOutName),'w')
        idsOutName = basename + ".ids"
        idsOutFile = file(os.path.join(outputLinkTrackletsInfilesDir, idsOutName),'w')
        writeDetsIdsFiles(detsOutFile, idsOutFile, allTrackletsFiles, allDias, cursor, diasDb, diasTable)
        detsOutFile.close()
        idsOutFile.close()
        # write C++ style start_t_range, which takes a fixed MJD.
        startTRangeOutCpp = os.path.join(outputLinkTrackletsInfilesDir, basename + ".date.start_t_range")
        startTRangeOutFileCpp = file(startTRangeOutCpp,'w')
        startTRangeOutFileCpp.write("%f"%(obsHistToExpMjd[obsHistsThisDataSet[-1]] + EPSILON))
        startTRangeOutFileCpp.close()
        

    

    if (WRITE_ADDITIONAL_METADATA):
        statsFilename = os.path.join(outputLinkTrackletsInfilesDir, basename + ".info")
        statsFile = file(statsFilename,'w')
        statsFile.write("!num_start_images num_support_images start_image_first_date start_image_last_date support_image_first_date support_image_last_date\n")
        statsFile.write("%d %d %f %f %f %f\n" % (len(obsHistsThisDataSet), len(supportObsHists),
                                                 obsHistToExpMjd[obsHistsThisDataSet[0]],
                                                 obsHistToExpMjd[obsHistsThisDataSet[-1]],
                                                 obsHistToExpMjd[supportObsHists[0]],
                                                 obsHistToExpMjd[supportObsHists[-1]]))

        for [obsHistCollection, headerString] in [[obsHistsThisDataSet, "FIRST ENDPOINT IMAGES"],
                                                  [supportObsHists, "SUPPORT IMAGES"]]:            
            statsFile.write("\n%s: by obsHistId fieldId expMjd trackletRootedInImage\n" % headerString)
            for obsHist in obsHistCollection:
                statsFile.write("%d %d %2.10f %d\n" % (obsHist, obsHistToFieldId[obsHist], obsHistToExpMjd[obsHist], obsHistToNumTlets[obsHist]))
        statsFile.close()

    print "finished writing output file for ", basename, " at ", time.ctime()
    sys.stdout.flush()
    

def appendSlashIfNeeded(dirName):                        
    if dirName[-1] != "/":
        dirName = dirName + "/"
    return dirName

if __name__=="__main__":

    trackletsByObsHistDir= appendSlashIfNeeded(sys.argv[1])

    # place to put .miti files for input to c linkTracklets
    outputLinkTrackletsInfilesDir=appendSlashIfNeeded(sys.argv[2])
    diasDb = sys.argv[3]
    diasTableOrDiasFile = sys.argv[4]
    trackingWindowDays = int(sys.argv[5])

    print "Reading tracklets from ", trackletsByObsHistDir
    print "Writing linkTracklets infiles to ", outputLinkTrackletsInfilesDir
    if PRELOAD_DIAS_FROM_FILE:
        print "Reading diaSources from: ", diasTableOrDiasFile
    else:
        print "Reading diaSources from: ", diasDb , ".", diasTable
    print "Reading image info from: ", mopsDatabases.OPSIM_DB , ".", mopsDatabases.OPSIM_TABLE
    print "Writing C-style MITI files: ", WRITE_C_INFILES
    print "Writing C++ style dets/ids files: ", WRITE_CPP_INFILES
    print "Preloading DiaSources into memory: ", PRELOAD_DIAS
    print "Number of days in tracking window: ", trackingWindowDays
    sys.stdout.flush()
    # SSCursor works much better for massive reads. It has a really ugly interface, though, which
    # doesn't work with the non-PRELOAD_DIAS code.
    if PRELOAD_DIAS:
        dbcurs = mopsDatabases.getCursor(useSSCursor=True)
    else:
        dbcurs = mopsDatabases.getCursor(useSSCursor=False)        
    makeLinkTrackletsInfiles(dbcurs, trackletsByObsHistDir, outputLinkTrackletsInfilesDir, 
                             diasDb, diasTableOrDiasFile, trackingWindowDays)
    print "DONE."

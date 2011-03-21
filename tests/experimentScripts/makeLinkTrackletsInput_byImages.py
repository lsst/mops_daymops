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

import mopsDatabases

# YOU ABSOLUTELY MUST SET THESE ARGUMENTS

TRACKLETS_BY_OBSHIST_DIR="/workspace1/jmyers/nightlyDiasAstromErr/tracklets/collapsed/byObsHistId/"

# place to put .miti files for input to c linkTracklets
OUTPUT_LINKTRACKLETS_INFILE_DIR="/workspace0/jmyers/nightlyDiasAstromErr_linkTrackletsInfiles_maxv0.5_mint_0.01_15dayWindows_CandCPP/"


# place to put start_t_ranges for each linkTracklets input files
OUTPUT_START_T_RANGE_FILES_DIR=OUTPUT_LINKTRACKLETS_INFILE_DIR



# obscode added to MITI files.
FORCED_OBSCODE="807"

# this needs to be larger than the min time between images.
EPSILON=1e-5

TRACKING_WINDOW_DAYS=15 # in days; we only look for tracks spanning <= this number of nights.
MAX_START_IMAGES_PER_RUN=30

TRACKLETS_BY_OBSHIST_SUFFIX=".tracklets.byDiaId"
OBSHIST_TO_TRACKLETS_FILE=lambda x: os.path.join(TRACKLETS_BY_OBSHIST_DIR) + str(x) + TRACKLETS_BY_OBSHIST_SUFFIX
TRACKLETS_FILE_TO_OBSHIST=lambda x: (int(os.path.basename(x)[:-len(TRACKLETS_BY_OBSHIST_SUFFIX)]))

# if true, will put some metadata about each data set into the start t range files directory.
WRITE_ADDITIONAL_METADATA=True

# jmyers: this seems to work for our data set but test your own first!
MJD_TO_NIGHT_NUM=lambda mjd: int(mjd)


# set to True if dias will fit in memory; they will load at beginning
# of execution and this will allow for MUCH faster lookups later.
# set to False if debugging to avoid a big wait at startup.
PRELOAD_DIAS=True




WRITE_CPP_INFILES=True

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




def fetchAllDiasFromDb(cursor):
    """ return a dictionary mapping diaId to diaSource, for all Dias
    in our data set."""
    toRet = {}
    print "Fetching DiaSources from DB into memory..."
    s = """ SELECT diaSourceId, ra, decl, taiMidPoint, mag, ssmId FROM
            %s.%s; """ % (mopsDatabases.DIAS_DB, mopsDatabases.DIAS_TABLE)
    cursor.execute(s)
    results = cursor.fetchall()
    for row in results:
        [diaId, ra, dec, mjd, mag, objId] = row
        d = Detection(diaId=diaId, ra=ra, dec=dec, mjd=mjd, mag=mag, objId=objId)
        toRet[diaId] = d

    print "... Done fetching dias."
    return toRet




def lookUpDia(allDias, diaId, cursor=None):
    if PRELOAD_DIAS:
        return allDias[diaId]
    else:
        s = """ SELECT diaSourceId, ra, decl, taiMidPoint, mag, ssmId FROM 
                %s.%s WHERE diaSourceId=%d; """ % (mopsDatabases.DIAS_DB, mopsDatabases.DIAS_TABLE, diaId)
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
    


def makeLinkTrackletsInfiles(cursor):
    # fetch all dias into memory (just use DB). This will be fasted in the long run than doing lots of    
    # cross-language/client-server calls to the DB when we need them.
    if PRELOAD_DIAS:
        allDias = fetchAllDiasFromDb(cursor)
    else:
        allDias = None

    # fetch all images which generated tracklets into memory (use a glob of the obshist_to_tracklets directory)
    # get MJDs of each image (use DB)
    trackletsFiles = glob.glob(os.path.join(TRACKLETS_BY_OBSHIST_DIR, '*' + TRACKLETS_BY_OBSHIST_SUFFIX))
    allObsHists = map(TRACKLETS_FILE_TO_OBSHIST, trackletsFiles)
    obsHistToExpMjd, obsHistToFieldId = getExpMjdAndFieldIdForObsHists(cursor, allObsHists)
                                     
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
        # now get those from nights n + 1 throught n + TRACKING_WINDOW_DAYS
        supportNights = range(curNightNum + 1, curNightNum+TRACKING_WINDOW_DAYS)
        for supportNight in supportNights:
            if nightNumsToObsHists.has_key(supportNight):
                supportObsHists += nightNumsToObsHists[supportNight]

        # if we have ANY support images at all, then go ahead and write output
        if len(supportObsHists) > 0:
            writeOutputFiles(allDias, cursor, obsHistsThisDataSet, supportObsHists, obsHistToExpMjd, obsHistToFieldId)




def writeMitiTrackletsToOutFile(mitiOut, mitiCheatSheetFile, allDias, cursor, trackletsFile, firstTrackletId):

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

            det = lookUpDia(allDias, diaId, cursor=cursor)
            mitiOut.write("%d %2.10f %2.10f %2.10f %2.10f %s %s 0.0 0.0\n" % \
                              (firstTrackletId + tletsWritten,
                               det.mjd, det.ra, det.dec, det.mag, FORCED_OBSCODE, 
                               det.objId))
            mitiCheatSheetFile.write("%d\n" % diaId)

        line = trackletsFile.readline()

        tletsWritten += 1

    return tletsWritten

    
def writeDetsIdsFiles(detsOutFile, idsOutFile, allTrackletsFileNames, allDias, cursor):
    """ detsOut and idsOut are expected to be open files where output is written.

    writes the simpler linkTracklets input as would be used for C++
    linkTracklets; dets are all the diaSources (miti format) for the
    data set and ids are the tracklets (as sets of LINE NUMBERS into
    the dtes file) for all tracklets in the data set."""

    # first figure out what diaIds need to be written.
    allIds = []
    for trackletsFileName in allTrackletsFileNames:
        trackletsFile = file(trackletsFileName, 'r')
        tletLine = trackletsFile.readline()
        while tletLine != "":
            ids = map(int, tletLine.split())
            allIds += ids
            tletLine = trackletsFile.readline()
        trackletsFile.close()

    # write dias that are needed.
    lineNum = 0
    diaToLineNum = {}
    for diaId in allIds:            
        det = lookUpDia(allDias, diaId, cursor=cursor)
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




def writeOutputFiles(allDias, cursor, obsHistsThisDataSet, supportObsHists, obsHistToExpMjd, obsHistToFieldId):

    """ write tracklets from the following obsHists into outfiles
    with locations specified by constants in this file."""    

    basename = "image_%d_through_%d" % (obsHistsThisDataSet[0], obsHistsThisDataSet[-1])

    print "Writing output file for ", basename

    mitiOutName = basename + ".miti"
    mitiCheatSheetName = basename + ".miti.diaIds"
    mitiOut = file(os.path.join(OUTPUT_LINKTRACKLETS_INFILE_DIR, mitiOutName),'w')
    mitiCheatSheetFile = file(os.path.join(OUTPUT_LINKTRACKLETS_INFILE_DIR, mitiCheatSheetName),'w')

    curTrackletId = 0
    obsHistToNumTlets = {}
    for obsHist in obsHistsThisDataSet + supportObsHists:
        trackletsFileName = OBSHIST_TO_TRACKLETS_FILE(obsHist)
        trackletsFile = file(trackletsFileName, 'r')
        tletsWritten = writeMitiTrackletsToOutFile(mitiOut, mitiCheatSheetFile, allDias, cursor, trackletsFile, curTrackletId)                
        curTrackletId += tletsWritten
        trackletsFile.close()
        obsHistToNumTlets[obsHist] = tletsWritten

    mitiOut.close()
    mitiCheatSheetFile.close()

    # turns out C linkTracklets takes start_t_range as a time offset (in days) from the first image. 
    startTRangeOut = os.path.join(OUTPUT_START_T_RANGE_FILES_DIR, basename + ".start_t_range")
    startTRangeOutFile = file(startTRangeOut,'w')
    startTRangeOutFile.write("%f"%(obsHistToExpMjd[obsHistsThisDataSet[-1]] - obsHistToExpMjd[obsHistsThisDataSet[0]] + EPSILON))
    startTRangeOutFile.close()
    
    

    if WRITE_CPP_INFILES:
        # new: write C++ style outputs as well.
        allTrackletsFiles = []
        for obsHist in obsHistsThisDataSet + supportObsHists:
            trackletsFileName = OBSHIST_TO_TRACKLETS_FILE(obsHist)
            allTrackletsFiles.append(trackletsFileName)

        detsOutName = basename + ".dets"
        detsOutFile = file(os.path.join(OUTPUT_LINKTRACKLETS_INFILE_DIR, detsOutName),'w')
        idsOutName = basename + ".ids"
        idsOutFile = file(os.path.join(OUTPUT_LINKTRACKLETS_INFILE_DIR, idsOutName),'w')
        writeDetsIdsFiles(detsOutFile, idsOutFile, allTrackletsFiles, allDias, cursor)
        detsOutFile.close()
        idsOutFile.close()
        # write C++ style start_t_range, which takes a fixed MJD.
        startTRangeOutCpp = os.path.join(OUTPUT_START_T_RANGE_FILES_DIR, basename + ".date.start_t_range")
        startTRangeOutFileCpp = file(startTRangeOutCpp,'w')
        startTRangeOutFileCpp.write("%f"%(obsHistToExpMjd[obsHistsThisDataSet[-1]] + EPSILON))
        startTRangeOutFileCpp.close()
        

    

    if (WRITE_ADDITIONAL_METADATA):
        statsFilename = os.path.join(OUTPUT_START_T_RANGE_FILES_DIR, basename + ".info")
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

    print "finished writing output file for ", basename
                                                    


if __name__=="__main__":
    dbcurs = mopsDatabases.getCursor()
    makeLinkTrackletsInfiles(dbcurs)

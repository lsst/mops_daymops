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
DIAS_SUFFIX=".dias"

TRACKLETS_FILE_TO_OBSHIST=lambda x: (int(os.path.basename(x)[:-len(TRACKLETS_BY_OBSHIST_SUFFIX)]))

# if true, will put some metadata about each data set into the start t range files directory.
WRITE_ADDITIONAL_METADATA=True

# jmyers: this seems to work for our data set but test your own first!
MJD_TO_NIGHT_NUM=lambda mjd: int(mjd)







def obsHistToTrackletsFile(obsHist, obsHistDir):
    return os.path.join(obsHistDir, str(obsHist) + TRACKLETS_BY_OBSHIST_SUFFIX)

def obsHistToDiasFile(obsHist, diasObsHistDir):
    return os.path.join(diasObsHistDir, str(obsHist) + DIAS_SUFFIX)




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



def makeLinkTrackletsInfiles(dbcurs, trackletsByObsHistDir, diasByObsHistDir,
                             outputLinkTrackletsInfilesDir, trackingWindowDays):

    # figure out IDs of all images which generated tracklets (use a
    # glob of the obshist_to_tracklets directory) and their MJDs (use
    # DB)
    trackletsFiles = glob.glob(os.path.join(trackletsByObsHistDir, 
                                            '*' + TRACKLETS_BY_OBSHIST_SUFFIX))
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
            writeOutputFiles(dbcurs, obsHistsThisDataSet, supportObsHists, obsHistToExpMjd,
                             obsHistToFieldId, nightNumsToObsHists,
                             outputLinkTrackletsInfilesDir, 
                             trackletsByObsHistDir, diasByObsHistDir)






def writeDetsIdsFiles(detsOutFile, idsOutFile, allTrackletsFileNames, 
                      diasByObsHistDir, neededImgDiaFiles):
    """ detsOut and idsOut are expected to be open files where output
    is to be written.

    writes the simpler linkTracklets input as would be used for C++
    linkTracklets; dets are all the diaSources for the
    data set and ids are the tracklets (as sets of LINE NUMBERS into
    the dets file) for all tracklets in the data set."""

    # write dets files while also tracking ids-> line numbers in our new file.
    curLine = 0
    diaToLineNum= {}
    for diasFile in neededImgDiaFiles:
        try:
            f = file(diasFile,'r')
            line = f.readline()
            while line != "":
                diaId = int(line.split()[0])
                diaToLineNum[diaId] = curLine
                detsOutFile.write(line)
                line = f.readline()
                curLine += 1
        except:
            # some images don't have associated diasources; this is
            # probably not a problem (especially when we are doing
            # simulations which only look at the ecliptic.)
            pass
    detsOutFile.close()

    # now translate from diaIds to indices (line nums) into the file.
    for trackletsFileName in allTrackletsFileNames:
        trackletsFile = file(trackletsFileName, 'r')
        tletLine = trackletsFile.readline()
        while tletLine != "":
            ids = map(int, tletLine.split())
            lineNums = []
            for diaId in ids:
                lineNums.append(diaToLineNum[diaId])
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


def getNeededImages(trackletObsHists, obsHistToExpMjd, cursor):
    """ figure out the earliest and latest obsHists in
    trackletObsHists.  Return every intermediate obsHist as well as
    any others from the last night."""
    tletMjds = [ obsHistToExpMjd[oh] for oh in trackletObsHists]
    earliestMjd = min(tletMjds)
    latestMjd = max(tletMjds)

    cursor.execute("""SELECT obsHistId FROM %s.%s WHERE expMJD >= %d -.001 AND 
expMJD < %d + .5;""" % (mopsDatabases.OPSIM_DB, mopsDatabases.OPSIM_TABLE, earliestMjd, latestMjd)) 
    res = cursor.fetchall()
    # we get e.g. [[obs1],[obs2]...] so flatten that
    realRes = map(lambda x: x[0], res)
    return realRes


def writeOutputFiles(cursor, obsHistsThisDataSet, supportObsHists, 
                     obsHistToExpMjd, obsHistToFieldId, nightNumToObsHists,
                     outputLinkTrackletsInfilesDir, 
                     trackletsByObsHistDir, diasByObsHistDir):

    """ write tracklets from the following obsHists into outfiles
    with locations specified by constants in this file."""

    basename = "image_%d_through_%d" % (obsHistsThisDataSet[0], obsHistsThisDataSet[-1])

    print "Writing output file for ", basename, " at ", time.ctime()
    sys.stdout.flush()

    # we need to count the number of
    # tracklets per image for our stats file.

    obsHistToNumTlets = {}
    for obsHist in obsHistsThisDataSet + supportObsHists:
        trackletsFileName = obsHistToTrackletsFile(obsHist, trackletsByObsHistDir)
        trackletsFile = file(trackletsFileName, 'r')
        tlets = countLines(trackletsFile)
        trackletsFile.close()
        obsHistToNumTlets[obsHist] = tlets


    needeTrackletsFiles = []
    for obsHist in obsHistsThisDataSet + supportObsHists:
        trackletsFileName = obsHistToTrackletsFile(obsHist, trackletsByObsHistDir)
        needeTrackletsFiles.append(trackletsFileName)

    # figure out what images' detections may be referenced by these
    # tracklets.  This should start with the first image of the
    # temporally earliest tracklet in the data set, through the last
    # image in which any tracklet started, plus any later image from
    # that same night as the latest tracklet.
    neededImages = getNeededImages(obsHistsThisDataSet + supportObsHists, 
                                   obsHistToExpMjd, cursor)
    neededImgDiaFiles = map(lambda x : obsHistToDiasFile(x, diasByObsHistDir), neededImages)

    detsOutName = basename + ".dets"
    detsOutFile = file(os.path.join(outputLinkTrackletsInfilesDir, detsOutName),'w')
    idsOutName = basename + ".ids"
    idsOutFile = file(os.path.join(outputLinkTrackletsInfilesDir, idsOutName),'w')
    writeDetsIdsFiles(detsOutFile, idsOutFile, needeTrackletsFiles, 
                      diasByObsHistDir, neededImgDiaFiles)
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
    diasByObsHistDir= appendSlashIfNeeded(sys.argv[2])

    # place to put .miti files for input to c linkTracklets
    outputLinkTrackletsInfilesDir=appendSlashIfNeeded(sys.argv[3])
    trackingWindowDays = int(sys.argv[4])

    print "Expect dias (grouped by image) in directory: ", diasByObsHistDir
    print "Expect tracklets (grouped by image) in directory: ", trackletsByObsHistDir
    print "Writing linkTracklets infiles to ", outputLinkTrackletsInfilesDir
    print "Reading image info from: ", mopsDatabases.OPSIM_DB , ".", mopsDatabases.OPSIM_TABLE
    print "Number of days in tracking window: ", trackingWindowDays

    sys.stdout.flush()

    dbcurs = mopsDatabases.getCursor(useSSCursor=False)
    makeLinkTrackletsInfiles(dbcurs, trackletsByObsHistDir, diasByObsHistDir,
                             outputLinkTrackletsInfilesDir, trackingWindowDays)
    print "DONE."

#!/usr/bin/env python7

"""

jmyers oct 26

Rather than build one linkTracklets input file per night, or make many
linkTracklets input files (one per N images) build linkTracklets input
files for one image.  This can be used to divide up the workload
between multiple machines on a cluster (like Abe).

Since there is no MySQL on Abe, this (optionally) version reads DiaSources/opSim
fields from a dump to local disk.  See "DIAS_FROM_FILE" below.

If reading a diaSource file, it should be in the format of Lynne's ephemeris dumps like 

MOPS64:/workspace1/jmyers/nightlyDiasAstromErr/dias_pt1_nodeep.short.astromErr.plusIds

The opsim dump should be like the output of dumpOpSimToFile.sql - it should be 

expMjd1 fieldId1 obsHistId1
expMjd2 fieldId2 obsHistId2
... etc.


"""

import os.path
import MySQLdb as db
import sys
import glob

# obscode added to MITI files.
FORCED_OBSCODE="807"

# this needs to be larger than the min time between images.
EPSILON=1e-5

TRACKING_WINDOW_DAYS=15 # in days; we only look for tracks spanning <= this number of nights.


# how to find tracklets (by Dia Id) organized by their "root" obsHistId.
TRACKLETS_BY_OBSHIST_DIR="/workspace1/jmyers/nightlyDiasAstromErr/tracklets/collapsed/byObsHistId/"
TRACKLETS_BY_OBSHIST_SUFFIX=".tracklets.byDiaId"
#these are not the dias in a given obsHist per se but the set of dias
#which are referenced by the corresponding tracklets
#file. buildDiaSourceSetsForTrackletFiles.py can create these files.
DIAS_BY_OBSHIST_DIR="/workspace1/jmyers/nightlyDiasAstromErr/tracklets/collapsed/byObsHistId/"
DIAS_BY_OBSHIST_SUFFIX=".tracklets.dias_lynneFormat"


OBSHIST_TO_TRACKLETS_FILE=lambda x: os.path.join(TRACKLETS_BY_OBSHIST_DIR) + str(x) + TRACKLETS_BY_OBSHIST_SUFFIX
TRACKLETS_FILE_TO_OBSHIST=lambda x: (int(os.path.basename(x)[:-len(TRACKLETS_BY_OBSHIST_SUFFIX)]))

OBSHIST_TO_DIAS_FILE=lambda x: os.path.join(DIAS_BY_OBSHIST_DIR) + str(x) + DIAS_BY_OBSHIST_SUFFIX
DIAS_FILE_TO_OBSHIST=lambda x: (int(os.path.basename(x)[:-len(DIAS_BY_OBSHIST_SUFFIX)]))

# if true, will put some metadata about each data set into the start t range files directory.
WRITE_ADDITIONAL_METADATA=True
#if true, write .ids and .dets files for C++ linkTracklets
WRITE_CPP_STYLE_INPUTS=False
#if true, write C-style MITI input
WRITE_C_STYLE_INPUTS=True

# jmyers: this seems to work for our data set but test your own first
MJD_TO_NIGHT_NUM=lambda mjd: int(mjd)


# set to True if dias will fit in memory and you want them read from
# the disk; they will load at beginning of execution and this will
# allow for faster lookups later.  set to False if debugging to avoid
# waiting on that all the time.  if set to False, it will read dias
# from MySQL and the DB stuff must be set below.
DIAS_FROM_FILE=True

DB_USER="jmyers"
DB_PASSWD="jmyers"
DB_HOST="localhost"
OPSIM_DB="opsim_3_61"
OPSIM_TABLE="output_opsim3_61"

#DIAS_DB="mops_noDeepAstromError"
#DIAS_TABLE="fullerDiaSource"





def readDias(diasDataFile, previouslyRead=None):

    """ reads a dump of dias, which include diaId expMjd ssmId
    obsHistId for every diaSource.  Returns a dict mapping diaId to
    data."""

    if previouslyRead == None:
        idToDias = {}
    else:
        idToDias = previouslyRead

    line = diasDataFile.readline()
    while line != "":
        [diaId, obsHistId, ssmId, ra, dec, expMjd, mag, snr] = line.split()
        diaId = int(diaId)
        [ra, dec, expMjd, mag, snr] = map(float, [ra, dec, expMjd, mag, snr])
        obsHistId = int(obsHistId)        
        
        idToDias[diaId] = DiaSource(diaId=diaId, obsTime=expMjd, ssmId=ssmId, obsHistId=obsHistId, 
                                    ra=ra, dec=dec, mag=mag)
        
        line = diasDataFile.readline()

    return idToDias





class DiaSource:
    def __init__(self, diaId, obsTime, ssmId, obsHistId, ra, dec, mag):
        self.diaId = diaId
        self.obsTime = obsTime
        self.ssmId = ssmId
        self.obsHistId = obsHistId
        self.ra = ra
        self.dec = dec
        self.mag = mag

    def getDiaId(self):
        return self.diaId

    def getObsTime(self):
        return self.obsTime

    def getSsmId(self):
        return self.ssmId

    def getObsHistId(self):
        return self.obsHistId
    
    def getRa(self):
        return self.ra

    def getDec(self):
        return self.dec

    def getMag(self):
        return self.mag




def lookUpDia(allDias, diaId, cursor=None):
    if DIAS_FROM_FILE:
        return allDias[diaId]
    else:
        s = """ SELECT diaSourceId, ra, decl, taiMidPoint, mag, ssmId FROM 
                %s.%s WHERE diaSourceId=%d; """ % (DIAS_DB, DIAS_TABLE, diaId)
        cursor.execute(s)
        [diaId, ra, dec, mjd, mag, objId] = cursor.fetchone()
        d = DiaSource(diaId=diaId, ra=ra, dec=dec, obsTime=mjd, mag=mag, objId=objId)
        return d



# TO DO: 

# need to implement below function call:

#  obsHistToExpMjd, obsHistToFieldId = getExpMjdAndFieldIdForObsHistsFromFile(opSimDumpFile, allObsHists)

def getExpMjdAndFieldIdForObsHistsFromFile(opSimDumpFileName, allObsHists):
    opSimDumpFile = file(opSimDumpFileName, 'r')
    line = opSimDumpFile.readline()
    obsHistToExpMjd = {}
    obsHistToFieldId = {}
    
    while line != "":
        [expMjd, fieldId, obsHistId] = line.split()
        expMjd = float(expMjd)
        [fieldId, obsHistId] = map(int, [fieldId, obsHistId])
        if obsHistId in allObsHists:
            obsHistToExpMjd[obsHistId] = expMjd
            obsHistToFieldId[obsHistId] = fieldId
        line = opSimDumpFile.readline()

    return obsHistToExpMjd, obsHistToFieldId


def getExpMjdAndFieldIdForObsHists(cursor, obsHists):

    """ use the DB to look up all the obsHists in obsHists and make a
    dictionary mapping obsHistId to time of image. Also make a dict
    mapping obsHistId to fieldId. Return both dictionaries in that
    order."""

    s = """ SELECT expMjd, fieldId, obsHistId FROM %s.%s where obsHistId in ( """ \
        % (OPSIM_DB, OPSIM_TABLE)

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



def getTrackletsFilesAndMetadata(trackletsInfilePattern, cursor=None, opSimDumpFile=None):

    if [cursor, opSimDumpFile] == [None, None]:
        raise Exception("makeLinkTrackletsInfile: must give either a DB connection for opsim or name of an input file containing the opsim records.")


    trackletsFiles = glob.glob(trackletsInfilePattern)
    allObsHists = map(TRACKLETS_FILE_TO_OBSHIST, trackletsFiles)

    if cursor != None:
        print "Reading ObsHist info from DB..."
        obsHistToExpMjd, obsHistToFieldId = getExpMjdAndFieldIdForObsHists(cursor, allObsHists)
        print "Done."
    else:
        print "Reading ObsHist info from file..."
        obsHistToExpMjd, obsHistToFieldId = getExpMjdAndFieldIdForObsHistsFromFile(opSimDumpFile, allObsHists)
        print "...Done."
           
    # bin images by night number
    allExpMjds = [ obsHistToExpMjd[oh] for oh in obsHistToExpMjd.keys() ]
    allNightNums = map(MJD_TO_NIGHT_NUM, allExpMjds)

    nightNumsToObsHists = {}
    for nn in allNightNums:
        nightNumsToObsHists[nn] = []
    for obsHist in obsHistToExpMjd.keys():
        nightNumsToObsHists[MJD_TO_NIGHT_NUM(obsHistToExpMjd[obsHist])].append(obsHist)
        
    return [trackletsFiles, allObsHists, nightNumsToObsHists, obsHistToExpMjd, obsHistToFieldId]



def getSupportImages(lastFirstEndpointObsHist, nightNumsToObsHists, obsHistToExpMjd):
    # DO THE ACTUAL WORK OF FINDING NEEDED OBSHISTS
    firstNightNum = MJD_TO_NIGHT_NUM(obsHistToExpMjd[lastFirstEndpointObsHist])
    supportObsHists = []
    
    for night in range(firstNightNum, firstNightNum + TRACKING_WINDOW_DAYS + 1):
        print "Looking for images on night ", night
        if nightNumsToObsHists.has_key(night):
            obsHistsTonight = nightNumsToObsHists[night]
            print "Found images with obsHistIds: ", obsHistsTonight
            for obsHist in obsHistsTonight:
                if obsHistToExpMjd[obsHist] > obsHistToExpMjd[lastFirstEndpointObsHist]:
                    print "  Adding ", obsHist, " at time ", obsHistToExpMjd[obsHist], "  to support set."
                    supportObsHists.append(obsHist)
    
    return supportObsHists


def makeLinkTrackletsInFile(startEndpointObsHist, outFileBaseName, cursor=None, diasInFile=None, opSimDumpFile=None):

    [trackletsFiles, allObsHists, nightNumToObsHists, obsHistToExpMjd, obsHistToFieldId] = \
        getTrackletsFilesAndMetadata(os.path.join(TRACKLETS_BY_OBSHIST_DIR, '*' + TRACKLETS_BY_OBSHIST_SUFFIX),
                                     cursor, opSimDumpFile)

    getSupportImages(startEndpointObsHist, nightNumToObsHists, obsHistToExpMjd)
    
    if len(supportObsHists) > 0:
        writeOutputFiles(outFileBaseName, cursor, [startEndpointObsHist], 
                         supportObsHists, obsHistToExpMjd, obsHistToFieldId)




def writeMitiTrackletsToOutFile(mitiOut, mitiCheatSheetFile, cursor, trackletsFile, diasFile, firstTrackletId):

    """ expects trackletsFile to be an open file full of tracklets in
    the diaId form (spaces between diaIds, newlines between
    tracklets).  mitiOut and mitiCheatSheetFile are expected to be an
    open file as well, where output is written.  curTrackletId is the
    first ID which will be used.  Returns the number of tracklets
    written."""

    tletsWritten = 0

    allDias = readDias(diasFile)

    line = trackletsFile.readline()
    while line != "":
        items = line.split()
        diaIds = map(int, items)
        for diaId in diaIds:
            #ID EPOCH_MJD RA_DEG DEC_DEG MAG OBSCODE OBJECT_NAME LENGTH ANGLE [EXPOSURE_TIME]

            det = lookUpDia(allDias, diaId, cursor=cursor)
            mitiOut.write("%d %2.10f %2.10f %2.10f %2.10f %s %s 0.0 0.0\n" % \
                              (firstTrackletId + tletsWritten,
                               det.getObsTime(), det.getRa(), det.getDec(), det.getMag(), FORCED_OBSCODE, 
                               det.getSsmId()))
            mitiCheatSheetFile.write("%d\n" % diaId)

        line = trackletsFile.readline()

        tletsWritten += 1

    return tletsWritten

    



def writeDetsIdsFiles(basename, allTrackletsAndDiasFileNames):
    """ 
    writes the simpler linkTracklets input as would be used for C++
    linkTracklets; dets are all the diaSources (miti format) for the
    data set and ids are the tracklets (as sets of diaIds, AND as sets
    of line numbers from the .dets file) for all tracklets in the data set.

    return a dict mapping obsHistId to number of tracklets rooted at the given obsHist.
    """
    obsHistToNumTlets = {}

    # first, figure out what detections are needed. Simultaneously
    # write the by-dia-id tracklets file.
    # also read the dets files and build the needed dict of DiaId -> detection
    # also count the number of tracklets at each image.
    
    idsOutFileByDiaId = file(basename + ".ids.byDiaId", 'w')
    allIds = set()
    allDias = {}
    for i in range(len(allTrackletsAndDiasFileNames)):

        trackletsFileName, diasFileName = allTrackletsAndDiasFileNames[i]
        diasFile = file(diasFileName,'r')
        allDias = readDias(diasFile, previouslyRead=allDias)
        diasFile.close()

        trackletsFile = file(trackletsFileName, 'r')
        tletLine = trackletsFile.readline()
        obsHistId = DIAS_FILE_TO_OBSHIST(diasFileName)
        numTlets = 0
        while tletLine != "":
            ids = map(int, tletLine.split())
            for diaId in ids:
                allIds.add(diaId)
            idsOutFileByDiaId.write(tletLine)
            tletLine = trackletsFile.readline()
            numTlets += 1
        trackletsFile.close()

        obsHistToNumTlets[obsHistId] = numTlets

    idsOutFileByDiaId.close()
    

    # write out all the needed dets, keep a map from ID to line number.
    detsOutFile = file(basename + ".dets",'w')
    detToLineNum = {}
    curLine = 0
    for diaId in allIds:
        det = lookUpDia(allDias, diaId)
        detsOutFile.write("%d %2.10f %2.10f %2.10f %2.10f %s %s 0.0 0.0\n" % \
                              (diaId,
                               det.getObsTime(), det.getRa(), det.getDec(), det.getMag(), FORCED_OBSCODE, 
                               det.getSsmId()))

        detToLineNum[diaId] = curLine
        curLine += 1
    
    # now we know the line numbers associated with each detection
    # written to the outfile.  make another pass over the tracklets files;
    # this time write the by-line-num linkTracklets in file for older
    # versions of linkTracklets.

    idsByLineNums = file(basename + ".ids.byLineNum",'w')

    for i in range(len(allTrackletsAndDiasFileNames)):
        trackletsFileName, diasFileName = allTrackletsAndDiasFileNames[i]
        trackletsFile = file(trackletsFileName, 'r')
        tletLine = trackletsFile.readline()
        while tletLine != "":
            ids = map(int, tletLine.split())
            lineNums = map(lambda x: detToLineNum[x], ids)
            
            for lineNum in lineNums:
                idsByLineNums.write("%d " % lineNum)
            idsByLineNums.write("\n")

            tletLine = trackletsFile.readline()

        trackletsFile.close()
    idsByLineNums.close()

    return obsHistToNumTlets




def writeOutputFiles(outFileBaseName, cursor, obsHistsThisDataSet, supportObsHists, obsHistToExpMjd, obsHistToFieldId):

    """ write tracklets from the following obsHists into outfiles
    with locations specified by constants in this file."""    

    basename = outFileBaseName

    print "Writing output file for ", basename

    if WRITE_C_STYLE_INPUTS:
        mitiOutName = basename + ".miti"
        mitiCheatSheetName = basename + ".miti.diaIds"
        mitiOut = file(mitiOutName,'w')
        mitiCheatSheetFile = file(mitiCheatSheetName,'w')

        curTrackletId = 0
        obsHistToNumTlets = {}
        for obsHist in obsHistsThisDataSet + supportObsHists:
            trackletsFileName = OBSHIST_TO_TRACKLETS_FILE(obsHist)
            trackletsFile = file(trackletsFileName, 'r')
            diasFileName = OBSHIST_TO_DIAS_FILE(obsHist)
            diasFile = file(diasFileName,'r')
            tletsWritten = writeMitiTrackletsToOutFile(mitiOut, mitiCheatSheetFile, cursor, trackletsFile, diasFile, curTrackletId)                
            curTrackletId += tletsWritten
            trackletsFile.close()
            diasFile.close()
            obsHistToNumTlets[obsHist] = tletsWritten

        mitiOut.close()
        mitiCheatSheetFile.close()

    # new: write C++ style outputs as well.
    if WRITE_CPP_STYLE_INPUTS:
        allTrackletsAndDiaFiles = []
        for obsHist in obsHistsThisDataSet + supportObsHists:
            trackletsFileName = OBSHIST_TO_TRACKLETS_FILE(obsHist)
            diasFileName = OBSHIST_TO_DIAS_FILE(obsHist)
            allTrackletsAndDiaFiles.append([trackletsFileName, diasFileName])

        detsOutName = basename + ".dets"
        detsOutFile = file(detsOutName,'w')
        idsOutName = basename + ".ids"
        idsOutFile = file(idsOutName,'w')
        obsHistToNumTlets = writeDetsIdsFiles(basename, allTrackletsAndDiaFiles)
        detsOutFile.close()
        idsOutFile.close()
    #end if WRITE_CPP_STYLE_INPUTS

    startTRangeOut =  basename + ".start_t_range"
    startTRangeOutFile = file(startTRangeOut,'w')
    startTRangeOutFile.write("%f"%(obsHistToExpMjd[obsHistsThisDataSet[-1]] + EPSILON))
    startTRangeOutFile.close()

    

    if (WRITE_ADDITIONAL_METADATA):
        statsFilename = os.path.join(basename + ".info")
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
    import sys

    if not DIAS_FROM_FILE:
        if len(sys.argv) != 3:
            raise Exception("USAGE: makeLinkTrackletsInput_perImage.py <obsHistId> <outFileName>")
        
        obsHistId = int(sys.argv[1])
        outFileName = sys.argv[2]

        dbconn = db.connect(host=DB_HOST, user=DB_USER, passwd=DB_PASSWD)
        dbcurs = dbconn.cursor()
        makeLinkTrackletsInfile(obsHistId, outFileName, cursor=dbcurs)

    else:
    
        if len(sys.argv) != 4:
            raise Exception("USAGE: makeLinkTrackletsInput_perImage.py <obsHistId> <outFileName> <opSimDump>")
        
        obsHistId = int(sys.argv[1])
        outFileName = sys.argv[2]
        opSimFileName = sys.argv[3]
        
        makeLinkTrackletsInFile(obsHistId, outFileName, opSimDumpFile=opSimFileName)

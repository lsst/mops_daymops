#!/usr/bin/env python


""" jmyers  nov 2 2010 

Some tools want input in the form of indices (line numbers) into a
file. Use this to get them.  I believe currently findTracklets writes
IDs but collapseTracklets and linkTracklets want indices into the
file.  Use this tool to convert from DiaIds to line numbers.

"""




import sys


def getIdToIndexMap(mitiInFile):
    line = mitiInFile.readline()
    idToIndex = {}
    lineNum = 0
    while line != "":
        mitiId = int(line.split()[0])
        idToIndex[mitiId] = lineNum
        lineNum += 1
        line = mitiInFile.readline()

    return idToIndex


def translateTracks(inFile, outFile, translationDict):
    line = inFile.readline()
    while line != "":
        items = map(int, line.split())
        for item in items:
            outFile.write("%d " % translationDict[item])
        outFile.write("\n")
        line = inFile.readline()


if __name__=="__main__":
    tracksIn = sys.argv[1]
    mitiIn = sys.argv[2]
    tracksOut = sys.argv[3]

    print "Detections are in ", mitiIn, ", reading linkages from ", tracksIn, " and writing to ", tracksOut


    tracksInFile = file(tracksIn, 'r')
    mitiInFile = file(mitiIn, 'r')
    tracksOutFile = file(tracksOut, 'w')

    idToIndex = getIdToIndexMap(mitiInFile)
    
    translateTracks(tracksInFile, tracksOutFile, idToIndex)
    tracksOutFile.close()
    print "Finished conversion successfully."

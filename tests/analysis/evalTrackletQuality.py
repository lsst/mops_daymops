#!/usr/bin/env python

""" jmyers 5/1/08

this script is for getting tracklet quality, average over all objects
and average by num of apparitions.

tracklet quality is:

 if tracklet is correct:  tracklet length / total object apparitions
 else:  0


where a non-lost apaprition is an apparition which appears in a
"correct" tracklet somewhere (where "correct" means "linking together
apparitions of exactly one object")
"""


import sys

from countObsPerObjects import trackletIsCorrect



def iterate(d, key):
    if d.has_key(key):
        d[key] += 1
    else:
        d[key] = 1



def trackletQuality(tracklet, assocIndices):
    return (1.0*len(tracklet))/(1.0*len(assocIndices))



def evalTrackletQuality(detections, pairs):
    objNameToTrackletsDict = {}
    objNameToIndicesDict = {}
    objNumObsToNumberDict = {}

    #find the indices associated with each object.
    
    for i in range(len(allDets)):
        objName = allDets[i][6]
        if objName != 'NS':
            if objNameToIndicesDict.has_key(objName):
                objNameToIndicesDict[objName].add(i)
            else:
                objNameToIndicesDict[objName] = set([i])

    #find number of objects which were seen X times for every possible X
    for name in objNameToIndicesDict.keys():
        iterate(objNumObsToNumberDict, len(objNameToIndicesDict[name]))




    sumQuality = 0.0
    objNumObsToSumQuality = {}
    objNumObsToNumTracklets = {}
    #objNumObsToNumTracklets[n] = number of *correct* tracklets associated with objects that had n apparitions

    for tracklet in allPairs:
        if (trackletIsCorrect(tracklet, detections)):
            objName = detections[tracklet[0]][6]
            assocIndices = objNameToIndicesDict[objName]
            tq = trackletQuality(tracklet, assocIndices)
            sumQuality += tq
            numApparitions = len(assocIndices)
            if objNumObsToSumQuality.has_key(numApparitions):
                objNumObsToSumQuality[numApparitions] += tq
                objNumObsToNumTracklets[numApparitions] += 1
            else:
                objNumObsToSumQuality[numApparitions] = tq
                objNumObsToNumTracklets[numApparitions] = 1
        #if tracklet is not correct, it has quality 0...
    
    averageQuality = 1.0*sumQuality / (1.0*len(pairs))
    aveQualityByNumObs = {}
    for key in objNumObsToNumberDict.keys():
        if objNumObsToSumQuality.has_key(key):
            aveQualityByNumObs[key] = float(objNumObsToSumQuality[key]) / \
                float(objNumObsToNumTracklets[key])
        else:
            aveQualityByNumObs[key] = 0.0

    return [averageQuality, aveQualityByNumObs]

if __name__=="__main__":
    detsFile = file(sys.argv[1], 'r')
    pairsFile = file(sys.argv[2], 'r')
    
    allDets = map(lambda x: x.split(), detsFile.readlines())
    allPairs = map(lambda x: map(lambda y: int(y), x.split()), pairsFile.readlines())

    
    [averageQuality, numApparitionsToQuality] = evalTrackletQuality(allDets, allPairs)

    print "Average Tracklet Quality: ", averageQuality
    for key in numApparitionsToQuality:
        print "Average Tracklet Quality, objects with ", key, " apparitions:\t", \
            numApparitionsToQuality[key]

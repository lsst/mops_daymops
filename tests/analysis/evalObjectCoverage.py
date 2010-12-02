#!/usr/bin/env python

""" jmyers 5/1/08

this script is for getting object coverage, average over all objects
and average by num of apparitions.

object coverage is:

 total unique non-lost apparitions of object / total object apparitions

where a non-lost apaprition is an apparition which appears in a
"correct" tracklet somewhere (where "correct" means "linking together
apparitions of exactly one object")
"""


import sys


from orderTracksByDet_byDiaId import Detection, readAllDets, lookUpDet



def trackletIsCorrect(trackletIds, allDets):
    objIds = set(map(lambda x: lookUpDet(allDets, x).objId, trackletIds))
    return (len(objIds) == 1 and (not '-1' in objIds) and (not "NS" in objIds))


def iterate(d, key):
    if d.has_key(key):
        d[key] += 1
    else:
        d[key] = 1



def objectCoverage(assocTracklets, assocIndices):
    uniqueIndices = set()
    for t in assocTracklets:
        for i in t:
            uniqueIndices.add(i)
    return (1.0*len(uniqueIndices))/(1.0*len(assocIndices))



def evalObjectCoverage(allDets, allTracklets):
    objNameToTrackletsDict = {}
    objNameToIdsDict = {}
    objNumObsToNumberDict = {}
    #find the indices associated with each object.
    
    for diaId in allDets.keys():
        objName = lookUpDet(allDets, diaId).objId
        if objName != 'NS' and objName != -1:
            if objNameToIdsDict.has_key(objName):
                objNameToIdsDict[objName].add(diaId)
            else:
                objNameToIdsDict[objName] = set([diaId])

    #find number of objects which were seen X times for every possible X
    for name in objNameToIdsDict.keys():
        iterate(objNumObsToNumberDict, len(objNameToIdsDict[name]))


    # associate correct tracklets with their respective objects 
    for tracklet in allTracklets:
        if trackletIsCorrect(tracklet, allDets):
            objName = lookUpDet(allDets, tracklet[0]).objId
            if objNameToTrackletsDict.has_key(objName):
                objNameToTrackletsDict[objName].append(tracklet)
            else:
                objNameToTrackletsDict[objName] = [ tracklet ]


    sumCoverage = 0.0
    objNumObsToNetCoverage = {}
    # find coverage percent of each object, add to sumCoverage,
    # add to objNumObsToNetCoverage

    for obj in objNameToTrackletsDict:
        objCo = objectCoverage(objNameToTrackletsDict[obj], objNameToIdsDict[obj])
        sumCoverage += objCo
        numObs = len(objNameToIdsDict[obj])
        if objNumObsToNetCoverage.has_key(numObs):
            objNumObsToNetCoverage[numObs] += objCo
        else:
            objNumObsToNetCoverage[numObs] = objCo

    # divide sumCoverage by total number of objects to get average coverage
    # divide each element of objNumObsToNetCoverage by the corresponding
    # objNumObsToNumberDict entry to get average coverage by object
    
    if len(objNameToTrackletsDict.keys()) > 0:
        averageCoverage = (1.0)*sumCoverage / (1.0*len(objNameToTrackletsDict.keys()))

        return [averageCoverage,len(objNameToIdsDict.keys())]

    else:
        return [0., 0]

if __name__=="__main__":
    detsFile = file(sys.argv[1], 'r')
    trackletsFile = file(sys.argv[2], 'r')
    
    allDets = readAllDets(detsFile)
    allTracklets = map(lambda x: map(lambda y: int(y), x.split()), trackletsFile.readlines())

    
    [averageCov, numObjects] = evalObjectCoverage(allDets, allTracklets)

    print "Average coverage: ", averageCov
    print "Num objects: ", numObjects
    #for key in numApparitionsToCov:
    #    print "Average coverage, objects with ", key, " apparitions:\t", \
    #        numApparitionsToCov[key]

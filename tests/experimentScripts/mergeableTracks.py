#!/usr/bin/env python

""" jmyers may 1

to help with debugging, find the tracks which could be merged into
longer tracks and report them.

Makes use of the output from readTracksWriteStats, which has the following format

[variable number of diaSources]  epoch:4.956344e+04 Best-fit RA p0, v, acc: 1.555300e+00 4.398283e-02 -1.671134e-01  Best-fit Dec p0, v, acc: -1.492982e+00 -5.426916e-02 -3.058764e-03 Underlying object: 1241985 Chi squared prob: Ra: 9.979549e-01 Dec: 9.998316e-01 fitRange: -4.660974e+00


"""



import sys

from binTrackletsByStartImage import getDiaTimesAndImages
import mopsDatabases

class TrackStats(object):
    def init(self, dias=None, epoch=None, ra0=None, raV=None, raAcc=None, 
             dec0=None, decV=None, decAcc=None, 
             objectName=None, probChisqRa=None, probChiqDec=None, fitRange=None):
        self.epoch = epoch
        self.ra0 = ra0
        self.raV = raV
        self.raAcc = raAcc
        self.dec0 = dec0
        self.decV = decV
        self.decAcc = decAcc
        self.objectName = objectName
        self.probChisqRa = probChisqRa
        self.probChisqDec = probChisqDec
        self.fitRange = fitRange
        self.dias = dias

    def fromString(self, s):
        try:
            items = s.split()
            self.dias = map(int, items[:-30])
            relevant = items[-29:]
            #print relevant
            self.epoch = float(relevant[0][6:])
            #print "epoch = ", epoch
            self.ra0, self.raV, self.raAcc = map(float, relevant[6:9])
            #print "Ra fit func = ",  self.ra0, self.raV, self.raAcc
            self.dec0, self.decV, self.decAcc = map(float, relevant[14:17])
            #print "dec func = ", self.dec0, self.decV, self.decAcc
            self.objectName = int(relevant[19])
            #print "name = ", self.objName
            self.probChisqRa = float(relevant[24])
            #print "probChisqR = ", self.probChisqRa
            self.probChisqDec = float(relevant[26])
            #print "probChisqD = ", self.probChisqDec
            self.fitRange = float(relevant[28])
            #print "fitRange = ", self.fitRange
            
        except:
            print "Couldn't parse string: ", s
            

def predPVatTime(track, predTime, raOrDec):
    if raOrDec == 'r':
        getP = lambda x: x.ra0
        getV = lambda x: x.raV
        getA = lambda x: x.raAcc
    elif raOrDec == 'd':
        getP = lambda x: x.dec0
        getV = lambda x: x.decV
        getA = lambda x: x.decAcc

    t0 = t.epoch
    v0 = getV(track)
    p0 = getP(track)
    acc = getA(track)

    dt = predTime - t.epoch
    vt = v0 + acc*(dt)
    pt = p0 + v0*dt + .2*acc*(dt**2)

    return pt, vt



def getUniqueNights(mjdsList):
    lastMjd = 0
    numNights = 0
    for mjd in sorted(mjdsList):
        if mjd - lastMjd > .5:
            numNights += 1
        lastMjd = mjd
    return numNights


if __name__=="__main__":
    inTracks = file(sys.argv[1])
    diasDb = sys.argv[2]
    diasTable = sys.argv[3]
    curs = mopsDatabases.getCursor()
    latestFirstEndpoint = float(sys.argv[4])
    print "Got latestFirstEndpoint = ", latestFirstEndpoint
    line = inTracks.readline()
    allDias = set()
    byObj = {}
    print "Reading in tracks..."
    while line != "":
        ts = TrackStats()
        ts.fromString(line)
        if ts.objectName != -1:
            if byObj.has_key(ts.objectName):
                byObj[ts.objectName].append(ts)
            else:
                byObj[ts.objectName] = [ts]
            allDias = allDias.union(set(ts.dias))
        line = inTracks.readline()

    numFindable = 0

    "Fetching ", len(allDias), " dias from DB..."
    diasTimesAndImages = getDiaTimesAndImages(allDias, curs,
                                              diasDb, diasTable)
    allTimes = [diasTimesAndImages[k][0] for k in diasTimesAndImages.keys()]
    firstTime = min(allTimes)
    lastTime = max(allTimes)
    aveTime = (firstTime + lastTime) / 2

    for obj in byObj:
        tracks = byObj[obj]
        epochs = set([t.epoch for t in tracks])
        
        if len(tracks) > 1:
            if min(epochs) < latestFirstEndpoint:
                objDias = set()
                for t in tracks:
                    objDias = objDias.union(set(t.dias))
                imgTimes = []
                for d in objDias:
                    time = diasTimesAndImages[d][0]
                    imgTimes.append(time)
                nights = getUniqueNights(list(imgTimes))
                otherNights = len(set(map(int, imgTimes)))
                          
                #if otherNights != nights:
                    #print " got nights = ", nights, " or ", otherNights 
                    #print "image times are ", imgTimes
                if (nights) >= 3:
                    numFindable += 1
                    raAccs = [t.raAcc for t in tracks]
                    decAccs = [t.decAcc for t in tracks]
                    print "max delta ra acc = ", (max(raAccs) - min(raAccs)) 
                    print "max delta dec acc = ", (max(decAccs) - min(decAccs)) 
                    # calculate predicted positions and velocities at aveTime
                    predPRas = []
                    predVRas = []
                    predPDecs = []
                    predVDecs = []
                    for track in tracks:
                        pr, vr = predPVatTime(track, aveTime, 'r')
                        pd, vd = predPVatTime(track, aveTime, 'd')
                        predPRas.append(pr)
                        predVRas.append(vr)
                        predPDecs.append(pd)
                        predVDecs.append(vd)

                    print "max delta ra V = ", (max(predVRas) - min(predVRas)) 
                    print "max delta dec V = ", (max(predVDecs) - min(predVDecs)) 
                    print "max delta ra pos = ", (max(predPRas) - min(predPRas)) 
                    print "max delta dec pos = ", (max(predPDecs) - min(predPDecs)) 
                    

                    #print "should get a mergeable track at  ",  min(epochs)
    print numFindable

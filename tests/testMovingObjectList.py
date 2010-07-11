#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import math
import random
import unittest

import MySQLdb as DBI

try:
    import lsst.mops.daymops.MovingObjectList as MovingObjectList
    from lsst.mops.daymops.MovingObject import MovingObject, STATUS
    from lsst.mops.daymops.Orbit import Orbit, STABLE_STATUS
    from lsst.mops.daymops.DiaSource import DiaSource
    from lsst.mops.daymops.Tracklet import Tracklet
except:
    raise(ImportError('Please setup daymops first.'))




# Constants.
DB_HOST = 'localhost'
DB_PORT = '3306'
DB_DBNAME = 'mops_onelunation'
DB_USER = 'jmyers'
DB_PASS = 'jmyers'




class TestMovingObjecList(unittest.TestCase):
    def checkObjs(self, obj1, obj2, attrName):
        fName = 'get%s%s' %(attrName[0].upper(), attrName[1:])
        v1 = getattr(obj1, fName)()
        v2 = getattr(obj2, fName)()
        
        if(isinstance(v1, float)):
            self.failUnlessAlmostEqual(v1, v2, 6,
                                       '%s is different: %s /= %s' %(attrName,
                                                                     str(v1),
                                                                     str(v2)))
        elif(isinstance(v1, list) or isinstance(v1, tuple)):
            # Lists are always floats in our case.
            for j in range(len(v1)):
                self.failUnlessAlmostEqual(v1[j], v2[j], 6,
                                       '%s is different: %s /= %s' %(attrName,
                                                                     str(v1),
                                                                     str(v2)))
        else:
            self.failUnlessEqual(v1, v2, '%s is different: %s /= %s' %(attrName,
                                                                       str(v1),
                                                                       str(v2)))
        return
        
        
    def setUp(self):
        # Build the connection string for the persistence framework.
        self.dbLocStr = 'mysql://%s:%d/%s' %(DB_HOST, int(DB_PORT), DB_DBNAME)
        
        # Connect to the database.
        self._dbh = DBI.connect(host=DB_HOST,
                                port=int(DB_PORT),
                                user=DB_USER,
                                passwd=DB_PASS,
                                db=DB_DBNAME)
        self._dbc = self._dbh.cursor()
        
        # Get all the Tracklets in the DB that are associated to MovingObjects.
        sql = '''\
select mops_Tracklet.trackletId, 
       mops_Tracklet.velRa, 
       mops_Tracklet.velDecl, 
       mops_Tracklet.velTot, 
       mops_Tracklet.status 
from mops_Tracklet
where mops_Tracklet.trackletId in 
(select mops_MovingObjectToTracklet.trackletId from mops_MovingObjectToTracklet)
order by mops_Tracklet.trackletId'''
        
        # Send the query.
        n = self._dbc.execute(sql)
        
        trueTracklets = {}
        row = self._dbc.fetchone()
        while(row):
            (_id, vRa, vDec, vTot, status) = row
            t = Tracklet()
            t.setTrackletId(int(_id))
            t.setVelRa(float(vRa))
            t.setVelDec(float(vDec))
            t.setVelTot(float(vTot))
            t.setStatus(str(status))
            trueTracklets[_id] = t
            
            row = self._dbc.fetchone()
        
        # Fetch trackletId -> diaSourceIds
        sql = '''\
select DIASource.diaSourceId, 
       DIASource.ra, 
       DIASource.decl, 
       DIASource.filterId, 
       DIASource.taiMidPoint, 
       DIASource.obsCode, 
       DIASource.apFlux, 
       DIASource.apFluxErr, 
       DIASource.refMag
from DIASource, mops_TrackletsToDIASource
where mops_TrackletsToDIASource.trackletId = %d and
      DIASource.diaSourceId = mops_TrackletsToDIASource.diaSourceId'''
        
        for _id in trueTracklets.keys():
            # Send the query.
            n = self._dbc.execute(sql %(_id))
            
            # Create the DiaSource objects manually.
            diaSources = []
            row = self._dbc.fetchone()
            while(row):
                (dId, ra, dec, fltr, mjd, ocode, fl, flErr, rMag) = row
                d = DiaSource()
                d.setDiaSourceId(int(dId))
                d.setRa(float(ra))
                d.setDec(float(dec))
                d.setFilterId(int(fltr))
                d.setTaiMidPoint(float(mjd))
                d.setObsCode(str(ocode))
                d.setApFlux(float(fl))
                d.setApFluxErr(float(flErr))
                d.setRefMag(float(rMag))
                diaSources.append(d)
                
                row = self._dbc.fetchone()
            
            # Add the diaSources to the Tracklet.
            trueTracklets[_id].setDiaSources(diaSources)
        
        # Now get the MovingObject -> Tracklet associations.
        trueMovingObjectToTracklets = {}
        sql = '''\
select movingObjectId, 
       trackletId 
from mops_MovingObjectToTracklet
order by movingObjectId'''
        n = self._dbc.execute(sql)
        if(not n):
                raise(DBSetupError('No MovingObject-Tracklet link!' %(_id)))
            
        row = self._dbc.fetchone()
        while(row):
            (moId, tId) = row
            if(trueMovingObjectToTracklets.has_key(moId)):
                trueMovingObjectToTracklets[moId].append(tId)
            else:
                trueMovingObjectToTracklets[moId] = [tId, ]
            row = self._dbc.fetchone()
        
        # Retrieve the list of MovingObject instances in the database.
        sql = '''\
select MovingObject.movingObjectId, 
       MovingObject.mopsStatus,
       MovingObject.h_v,
       MovingObject.g,
       MovingObject.arcLengthDays,
       MovingObject.q,
       MovingObject.e,
       MovingObject.i,
       MovingObject.node,
       MovingObject.argPeri,
       MovingObject.timePeri,
       MovingObject.epoch,
       MovingObject.src01,
       MovingObject.src02,
       MovingObject.src03,
       MovingObject.src04,
       MovingObject.src05,
       MovingObject.src06,
       MovingObject.src07,
       MovingObject.src08,
       MovingObject.src09,
       MovingObject.src10,
       MovingObject.src11,
       MovingObject.src12,
       MovingObject.src13,
       MovingObject.src14,
       MovingObject.src15,
       MovingObject.src16,
       MovingObject.src17,
       MovingObject.src18,
       MovingObject.src19,
       MovingObject.src20,
       MovingObject.src21,
       MovingObject.stablePass,
       MovingObject.meanAnom
from MovingObject'''
        
        # Send the query.
        n = self._dbc.execute(sql)
        
        # Create the DiaSource objects manually.
        self.trueMovingObjects = {}
        row = self._dbc.fetchone()
        while(row):
            src = [None, ] * 21
            (_id, 
             mopsStatus,
             h_v,
             g,
             arcLengthDays,
             q,
             e,
             i,
             node,
             argPeri,
             timePeri,
             epoch,
             src[0],
             src[1],
             src[2],
             src[3],
             src[4],
             src[5],
             src[6],
             src[7],
             src[8],
             src[9],
             src[10],
             src[11],
             src[12],
             src[13],
             src[14],
             src[15],
             src[16],
             src[17],
             src[18],
             src[19],
             src[20],
             stablePass,
             meanAnom) = row
            # Create the MovingObject and its Orbit instance.
            o = Orbit()
            o.setQ(float(q))
            o.setE(float(e))
            o.setI(float(i))
            o.setNode(node)
            o.setArgPeri(argPeri)
            o.setTimePeri(float(timePeri))
            o.setEpoch(float(epoch))
            o.setSrc(src)
            o.setStablePass(str(stablePass))
            o.setMeanAnom(float(meanAnom))
            
            m = MovingObject()
            m.setMovingObjectId(int(_id))
            m.setMopsStatus(str(mopsStatus))
            if(h_v != None):
                m.setH_v(float(h_v))
            else:
                m.setH_v(h_v)
            if(g != None):
                m.setG(float(g))
            else:
                m.setG(g)
            if(arcLengthDays != None):
                m.setArcLength(float(arcLengthDays))
            else:
                m.setArcLength(arcLengthDays)
            m.setOrbit(o)
            
            self.trueMovingObjects[_id] = m
            row = self._dbc.fetchone()
        
            # Manualy update the true MovingObjects. Only update those that
            # have Tracklets, clearly.
            if(_id in trueMovingObjectToTracklets.keys()):
                ts = [trueTracklets[tId] 
                      for tId in trueMovingObjectToTracklets[_id]]
                self.trueMovingObjects[_id].setTracklets(ts)
        return
    
    def testGetAllMovingObjectsShallow(self):
        # This actually means (see the function docs in MovingObjectList.py) 
        # retrieving all MovingObject with status /= merged.
        iter = MovingObjectList.getAllMovingObjects(self.dbLocStr,
                                                    shallow=True,
                                                    sliceId=0,
                                                    numSlices=1)
        movingObjects = dict([(m.getMovingObjectId(), m) for m in iter])
        
        trueNonMergedObjs = dict([(m.getMovingObjectId(), m) \
                                  for m in self.trueMovingObjects.values() \
                                  if m.getMopsStatus != STATUS['MERGED']])
        
        # Make sure that we got the same number of objects.
        self.failUnlessEqual(len(trueNonMergedObjs), len(movingObjects),
                             'number of retrieved objects is incorrect: %d /= %d' %(len(trueNonMergedObjs), len(movingObjects)))
        
        # Check 'em piece by piece. Do not test arcLength since we are not 
        # retrieving Tracklets.
        for _id in trueNonMergedObjs.keys():
            for attr in ('movingObjectId', 'mopsStatus', 'h_v', 'g'):
                self.checkObjs(trueNonMergedObjs[_id], movingObjects[_id], attr)
            for attr in ('q', 'e', 'i', 'node', 'argPeri', 'timePeri', 'epoch',
                         'meanAnom', 'src'):
                self.checkObjs(trueNonMergedObjs[_id].getOrbit(), 
                               movingObjects[_id].getOrbit(), attr)
        return
    
    def testGetAllMovingObjectsDeep(self):
        # This actually means (see the function docs in MovingObjectList.py) 
        # retrieving all MovingObject with status /= merged.
        iter = MovingObjectList.getAllMovingObjects(self.dbLocStr,
                                                    shallow=False,
                                                    sliceId=0,
                                                    numSlices=1)
        movingObjects = dict([(m.getMovingObjectId(), m) for m in iter])
        
        trueNonMergedObjs = dict([(m.getMovingObjectId(), m) \
                                  for m in self.trueMovingObjects.values() \
                                  if m.getMopsStatus != STATUS['MERGED']])
        
        # Make sure that we got the same number of objects.
        self.failUnlessEqual(len(trueNonMergedObjs), len(movingObjects),
                             'number of retrieved objects is incorrect: %d /= %d' %(len(trueNonMergedObjs), len(movingObjects)))
        
        # Check 'em piece by piece.
        for _id in trueNonMergedObjs.keys():
            for attr in ('movingObjectId', 'mopsStatus', 'h_v', 'g', 
                         'arcLength'):
                self.checkObjs(trueNonMergedObjs[_id], movingObjects[_id], attr)
            for attr in ('q', 'e', 'i', 'node', 'argPeri', 'timePeri', 'epoch',
                         'meanAnom', 'src'):
                self.checkObjs(trueNonMergedObjs[_id].getOrbit(), 
                               movingObjects[_id].getOrbit(), attr)
            
            # Check their Tracklets.
            trueTracklets = trueNonMergedObjs[_id].getTracklets()
            tracklets = movingObjects[_id].getTracklets()
            trueTracklets.sort()
            tracklets.sort()
            
            # Same numbers?
            self.failUnlessEqual(len(trueTracklets), len(tracklets),
                             'number of retrieved Tracklets for object %d is incorrect: %d /= %d' %(_id, len(trueTracklets), len(tracklets)))
            
            # Same Tracklets?
            for i in range(len(tracklets)):
                # Check Tracklet attributes.
                for attr in ('trackletId', 'velRa', 'velDec', 'status'):
                    self.checkObjs(trueTracklets[i], tracklets[i], attr)
        return
    
    def testGetAllUnstableMovingObjectsShallow(self):
        iter = MovingObjectList.getAllUnstableMovingObjects(self.dbLocStr,
                                                            shallow=True,
                                                            sliceId=0,
                                                            numSlices=1)
        movingObjects = dict([(m.getMovingObjectId(), m) for m in iter])
        
        trueNonMergedObjs = dict([(m.getMovingObjectId(), m) \
            for m in self.trueMovingObjects.values() \
            if m.getMopsStatus != STATUS['MERGED'] and \
               m.getOrbit().getStablePass() != STABLE_STATUS['STABLE']])
        
        # Make sure that we got the same number of objects.
        self.failUnlessEqual(len(trueNonMergedObjs), len(movingObjects),
                             'number of retrieved objects is incorrect: %d /= %d' %(len(trueNonMergedObjs), len(movingObjects)))
        
        # Check 'em piece by piece. Do not check arcLength since we are not 
        # retrieving Tracklets.
        for _id in trueNonMergedObjs.keys():
            for attr in ('movingObjectId', 'mopsStatus', 'h_v', 'g'):
                self.checkObjs(trueNonMergedObjs[_id], movingObjects[_id], attr)
            for attr in ('q', 'e', 'i', 'node', 'argPeri', 'timePeri', 'epoch',
                         'meanAnom', 'src'):
                self.checkObjs(trueNonMergedObjs[_id].getOrbit(), 
                               movingObjects[_id].getOrbit(), attr)
        return
    
    def testGetAllUnstableMovingObjectsDeep(self):
        iter = MovingObjectList.getAllUnstableMovingObjects(self.dbLocStr,
                                                            shallow=False,
                                                            sliceId=0,
                                                            numSlices=1)
        movingObjects = dict([(m.getMovingObjectId(), m) for m in iter])
        
        trueNonMergedObjs = dict([(m.getMovingObjectId(), m) \
            for m in self.trueMovingObjects.values() \
            if m.getMopsStatus != STATUS['MERGED'] and \
               m.getOrbit().getStablePass() != STABLE_STATUS['STABLE']])
        
        # Make sure that we got the same number of objects.
        self.failUnlessEqual(len(trueNonMergedObjs), len(movingObjects),
                             'number of retrieved objects is incorrect: %d /= %d' %(len(trueNonMergedObjs), len(movingObjects)))
        
        # Check 'em piece by piece. 
        for _id in trueNonMergedObjs.keys():
            for attr in ('movingObjectId', 'mopsStatus', 'h_v', 'g', 
                         'arcLength'):
                self.checkObjs(trueNonMergedObjs[_id], movingObjects[_id], attr)
            for attr in ('q', 'e', 'i', 'node', 'argPeri', 'timePeri', 'epoch',
                         'meanAnom', 'src'):
                self.checkObjs(trueNonMergedObjs[_id].getOrbit(), 
                               movingObjects[_id].getOrbit(), attr)
            
            # Check their Tracklets.
            trueTracklets = trueNonMergedObjs[_id].getTracklets()
            tracklets = movingObjects[_id].getTracklets()
            trueTracklets.sort()
            tracklets.sort()
            
            # Same numbers?
            self.failUnlessEqual(len(trueTracklets), len(tracklets),
                             'number of retrieved Tracklets for object %d is incorrect: %d /= %d' %(_id, len(trueTracklets), len(tracklets)))
            
            # Same Tracklets?
            for i in range(len(tracklets)):
                # Check Tracklet attributes.
                for attr in ('trackletId', 'velRa', 'velDec', 'status'):
                    self.checkObjs(trueTracklets[i], tracklets[i], attr)
        return






class DBSetupError(Exception): pass


def sanityCheck():
    # Make sure that the database and tables are there.
    dbh = None
    dbc = None
    try:
        dbh = DBI.connect(host=DB_HOST,
                          port=int(DB_PORT),
                          user=DB_USER,
                          passwd=DB_PASS,
                          db=DB_DBNAME)
        dbc = dbh.cursor()
    except:
        raise(DBSetupError('Please create and populate the %s database.' \
                           %(DB_DBNAME)))
    
    sql = 'select count(*) from mops_Tracklet'
    n = dbc.execute(sql)
    n = dbc.fetchone()[0]
    if(not n):
        raise(DBSetupError('Please populate the mops_Tracklet table.'))
    
    sql = 'select count(*) from MovingObject'
    n = dbc.execute(sql)
    n = dbc.fetchone()[0]
    if(not n):
        raise(DBSetupError('Please populate the MovingObject table.'))
    
    sql = 'select count(*) from mops_MovingObjectToTracklet'
    n = dbc.execute(sql)
    n = dbc.fetchone()[0]
    if(not n):
        raise(DBSetupError('Please populate the mops_MovingObjectToTracklet table.'))
    return



if(__name__ == '__main__'):
    sanityCheck()        
    unittest.main()





















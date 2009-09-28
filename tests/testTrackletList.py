#!/usr/bin/env python
import math
import random
import unittest

import MySQLdb as DBI

try:
    import lsst.mops.daymops.TrackletList as TrackletList
    from lsst.mops.daymops.Tracklet import Tracklet, STATUS
    from lsst.mops.daymops.DiaSource import DiaSource
except:
    raise(ImportError('Please setup daymops first.'))




# Constants.
DB_HOST = 'localhost'
DB_PORT = '3306'
DB_DBNAME = 'mops_test'
DB_USER = 'www'
DB_PASS = 'zxcvbnm'




class TestTrackletList(unittest.TestCase):
    def checkObjs(self, obj1, obj2, attrName):
        fName = 'get%s%s' %(attrName[0].upper(), attrName[1:])
        v1 = getattr(obj1, fName)()
        v2 = getattr(obj2, fName)()
        if(isinstance(v1, float)):
            self.failUnlessAlmostEqual(v1, v2, 6,
                                       '%s is different: %s /= %s' %(attrName,
                                                                     str(v1),
                                                                     str(v2)))
            return
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
        
        # Retrieve the list of DiaSource instances in the database as well as 
        # the ID of the DiaSources belonging to the current night.
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
from DIASource'''
        
        # Send the query.
        n = self._dbc.execute(sql)
        
        # Create the DiaSource objects manually.
        self.trueDiaSources = {}
        row = self._dbc.fetchone()
        while(row):
            (_id, ra, dec, fltr, mjd, ocode, fl, flErr, rMag) = row
            d = DiaSource()
            d.setDiaSourceId(int(_id))
            d.setRa(float(ra))
            d.setDec(float(dec))
            d.setFilterId(int(fltr))
            d.setTaiMidPoint(float(mjd))
            d.setObsCode(str(ocode))
            d.setApFlux(float(fl))
            d.setApFluxErr(float(flErr))
            d.setRefMag(float(rMag))
            self.trueDiaSources[_id] = d
            
            row = self._dbc.fetchone()
        
        
        # Then get all the Tracklets in the DB.
        # Retrieve the list of DiaSource instances in the database as well as 
        # the ID of the DiaSources belonging to the current night.
        sql = '''\
select mops_Tracklet.trackletId, 
       mops_Tracklet.velRa, 
       mops_Tracklet.velDecl, 
       mops_Tracklet.velTot, 
       mops_Tracklet.status 
from mops_Tracklet
order by mops_Tracklet.trackletId'''
        
        # Send the query.
        n = self._dbc.execute(sql)
        
        # Create the DiaSource objects manually.
        self.trueTracklets = {}
        row = self._dbc.fetchone()
        while(row):
            (_id, vRa, vDec, vTot, status) = row
            t = Tracklet()
            t.setTrackletId(int(_id))
            t.setVelRa(float(vRa))
            t.setVelDec(float(vDec))
            t.setVelTot(float(vTot))
            t.setStatus(str(status))
            self.trueTracklets[_id] = t
            
            row = self._dbc.fetchone()

        
        # Fetch trackletId -> diaSourceIds
        self.trueTrackletToDiaSources = {}
        sql = '''\
select diaSourceId 
from mops_TrackletsToDIASource 
where trackletId = %d'''
        for _id in self.trueTracklets.keys():
            n = self._dbc.execute(sql %(_id))
            if(not n):
                raise(DBSetupError('Tracklet %d has no DiaSources!' %(_id)))
            dIds = self._dbc.fetchall()
            self.trueTrackletToDiaSources[_id] = zip(*dIds)[0]
        
        
        # Now get the IDs of the ones for tonight.
        sql = 'select DIASourceIDTonight.DIASourceId from DIASourceIDTonight'
        n = self._dbc.execute(sql)
        self.trueDiaSourceIdsForTonight = [int(r[0]) for r in self._dbc.fetchall()]
        
        # Finally, get the MovingObject -> Tracklet associations.
        self.trueMovingObjectToTracklets = {}
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
            if(self.trueMovingObjectToTracklets.has_key(moId)):
                self.trueMovingObjectToTracklets[moId].append(tId)
            else:
                self.trueMovingObjectToTracklets[moId] = [tId, ]
            row = self._dbc.fetchone()
        return
    
    def testNewTrackletsFromTonight(self):
        # get the ids of tracklets that have at least one DiaSource ID in 
        # self.trueDiaSourceIdsForTonight
        tonight = set(self.trueDiaSourceIdsForTonight)
        tIds = [i for i in self.trueTrackletToDiaSources.keys() \
                if set(self.trueTrackletToDiaSources[i]).intersection(tonight)]
        
        # Now use TrackletList.newTrackletsFromTonight() and see what's up.
        iter = TrackletList.newTrackletsFromTonight(self.dbLocStr,
                                                    shallow=True,
                                                    sliceId=0,
                                                    numSlices=1)
        myTracklets = [t for t in iter]
        self.failUnlessEqual(len(tIds), len(myTracklets), 'number of tracklets retrieved is different: %d /= %d.' %(len(tIds), len(myTracklets)))
        
        # And now check them one by one. We only check Ids, since full retrieval
        # is tested elsewhere.
        for t in myTracklets:
            self.failUnless(t.getTrackletId() in tIds, 'Tracklet %d was incorrectly retrieved.' %(t.getTrackletId()))
        return
    
    def testMewTracklets(self):
        # Fetch all unattributed tracklets.
        tracklets = [t for t in TrackletList.newTracklets(self.dbLocStr,
                                                          shallow=False,
                                                          sliceId=0,
                                                          numSlices=1)]
        
        trueTracklets = [t for t in self.trueTracklets.values() \
                         if t.getStatus() == STATUS['UNATTRIBUTED']]
        self.failUnlessEqual(len(trueTracklets), len(tracklets), 'number of tracklets retrieved is different: %d /= %d.' %(len(trueTracklets), len(tracklets)))
        
        # Now check them one by one. It is a deep fetch, so check DiaSources as
        # well.
        for i in range(len(tracklets)):
            # Check Tracklet attributes.
            for attr in ('trackletId', 'velRa', 'velDec', 'status'):
                self.checkObjs(trueTracklets[i], tracklets[i], attr)
            
            # Now the DiaSources.
            trueSourcesIds = self.trueTrackletToDiaSources[trueTracklets[i].getTrackletId()]
            trueSources = [self.trueDiaSources[_id] for _id in trueSourcesIds]
            sources = tracklets[i].getDiaSources()
            
            self.failUnlessEqual(len(trueSources), len(sources), 'number of DiaSources is different: %d /= %d' %(len(trueSources), len(sources)))
            
            # And now the attributes.
            sources.sort()
            trueSources.sort()
            for j in range(len(sources)):
                for attr in ('diaSourceId', 'ra', 'dec', 'filterId', 
                             'taiMidPoint', 'obsCode', 'apFlux', 'apFluxErr', 
                             'refMag'):
                    self.checkObjs(trueSources[j], sources[j], attr)
        return
    
    def testTrackletsForMovingObject(self):
        # For each MovingObjectId, make sure that we get all the corresponding
        # Tracklets.
        for movingObjectId in self.trueMovingObjectToTracklets.keys():
            trueTrackletIds = self.trueMovingObjectToTracklets[movingObjectId]
            
            # Fetch the Tracklets.
            tIter = TrackletList.allTrackletsForMovingObject(self.dbLocStr,
                                                             movingObjectId,
                                                             True,
                                                             0,
                                                             1)
            trackletIds = [t.getTrackletId() for t in tIter]
            self.failUnlessEqual(len(trueTrackletIds), len(trackletIds), 'number of tracklets retrieved for MovingObject %d is different: %d /= %d.' %(movingObjectId, len(trueTrackletIds), len(trackletIds)))
            
            trueTrackletIds.sort()
            trackletIds.sort()
            for i in range(len(trueTrackletIds)):
                self.failUnlessEqual(trueTrackletIds[i], trackletIds[i], 'IDs of Tracklet number %d differ: %d /= %d' %(i, trueTrackletIds[i], trackletIds[i]))
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
    
    sql = 'select count(*) from DIASource'
    n = dbc.execute(sql)
    n = dbc.fetchone()[0]
    if(not n):
        raise(DBSetupError('Please populate the DIASource table.'))
    
    sql = 'select count(*) from DIASourceIDTonight'
    n = dbc.execute(sql)
    n = dbc.fetchone()[0]
    if(not n):
        raise(DBSetupError('Please populate the DIASourceIDTonight table.'))
    return



if(__name__ == '__main__'):
    sanityCheck()        
    unittest.main()





















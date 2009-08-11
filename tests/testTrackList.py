#!/usr/bin/env python
import math
import random
import unittest

import MySQLdb as DBI

try:
    import lsst.mops.daymops.TrackList as TrackList
    from lsst.mops.daymops.Track import Track
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




class TestTrackList(unittest.TestCase):
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
        
        # Fetch trackId -> trackletIds
        self.trueTrackToTracklets = {}
        sql = '''\
select trackId, trackletId 
from  mops_TracksToTracklet  
order by trackId'''
        n = self._dbc.execute(sql)
        if(not n):
            raise(DBSetupError('No Tracks in the database!'))
        row = self._dbc.fetchone()
        while(row):
            if(row[0] in self.trueTrackToTracklets.keys()):
                self.trueTrackToTracklets[row[0]].append(int(row[1]))
            else:
                self.trueTrackToTracklets[row[0]] = [int(row[1]), ]
            row = self._dbc.fetchone()
        return
    
    def testNewTracksShallow(self):
        iter = TrackList.newTracks(self.dbLocStr,
                                   shallow=True,
                                   sliceId=0,
                                   numSlices=1)
        tracks = dict([(t.getTrackId(), t) for t in iter])
        
        # Same number?
        self.failUnlessEqual(len(self.trueTrackToTracklets), len(tracks),
                             'the number of Tracks is different: %d /= %d' %(len(self.trueTrackToTracklets), len(tracks)))
        
        # And now, one by one.
        for trackId in self.trueTrackToTracklets.keys():
            # TrackId has to be there.
            self.failUnless(trackId in tracks.keys())
            
            track = tracks[trackId]
            trueTrack = self.trueTrackToTracklets[trackId]
            
            # And now check the Tracklets.
            trueTracklets = [self.trueTracklets[tId] for tId in trueTrack]
            tracklets = track.getTracklets()
            
            self.failUnlessEqual(len(trueTracklets), len(tracklets),
                                 'the number of Tracklets is different: %d /= %d' %(len(trueTracklets), len(tracklets)))
            if(not trueTracklets):
                continue
            
            # This is a deep fetch.
            for i in range(len(trueTracklets)):
                # Check Tracklet attributes.
                for attr in ('trackletId', ):
                    self.checkObjs(trueTracklets[i], tracklets[i], attr)
        return
    
    
    def testNewTracksDeep(self):
        iter = TrackList.newTracks(self.dbLocStr,
                                   shallow=False,
                                   sliceId=0,
                                   numSlices=1)
        tracks = dict([(t.getTrackId(), t) for t in iter])
        
        # Same number?
        self.failUnlessEqual(len(self.trueTrackToTracklets), len(tracks),
                             'the number of Tracks is different: %d /= %d' %(len(self.trueTrackToTracklets), len(tracks)))
        
        # And now, one by one.
        for trackId in self.trueTrackToTracklets.keys():
            # TrackId has to be there.
            self.failUnless(trackId in tracks.keys())
            
            track = tracks[trackId]
            trueTrack = self.trueTrackToTracklets[trackId]
            
            # And now check the Tracklets.
            trueTracklets = [self.trueTracklets[tId] for tId in trueTrack]
            tracklets = track.getTracklets()
            
            self.failUnlessEqual(len(trueTracklets), len(tracklets),
                                 'the number of Tracklets is different: %d /= %d' %(len(trueTracklets), len(tracklets)))
            if(not trueTracklets):
                continue
            
            # This is a deep fetch.
            for i in range(len(trueTracklets)):
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
                
    
    def testNewTracksDeepParallel(self):
        # Get half of the Tracks.
        iter = TrackList.newTracks(self.dbLocStr,
                                   shallow=False,
                                   sliceId=0,
                                   numSlices=2)
        tracks1 = dict([(t.getTrackId(), t) for t in iter])
        
        # Get the other half.
        iter = TrackList.newTracks(self.dbLocStr,
                                   shallow=False,
                                   sliceId=1,
                                   numSlices=2)
        tracks2 = dict([(t.getTrackId(), t) for t in iter])
        
        # Make sure that we do not get any duplicates.
        for _id in tracks1.keys():
            self.failUnless(_id not in tracks2.keys(),
                            'Track %d is present in both dictionaries!' %(_id))
        
        # Now concatenate the two dictionaries.
        tracks1.update(tracks2)
        tracks = tracks1
        
        # Same number?
        self.failUnlessEqual(len(self.trueTrackToTracklets), len(tracks),
                             'the number of Tracks is different: %d /= %d' %(len(self.trueTrackToTracklets), len(tracks)))
        
        # And now, one by one.
        for trackId in self.trueTrackToTracklets.keys():
            # TrackId has to be there.
            self.failUnless(trackId in tracks.keys())
            
            track = tracks[trackId]
            trueTrack = self.trueTrackToTracklets[trackId]
            
            # And now check the Tracklets.
            trueTracklets = [self.trueTracklets[tId] for tId in trueTrack]
            tracklets = track.getTracklets()
            
            self.failUnlessEqual(len(trueTracklets), len(tracklets),
                                 'the number of Tracklets is different: %d /= %d' %(len(trueTracklets), len(tracklets)))
            if(not trueTracklets):
                continue
            
            # This is a deep fetch.
            for i in range(len(trueTracklets)):
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
    
    sql = 'select count(*) from mops_TracksToTracklet'
    n = dbc.execute(sql)
    n = dbc.fetchone()[0]
    if(not n):
        raise(DBSetupError('Please populate the mops_TracksToTracklet table.'))
    return
    
    sql = 'select count(*) from mops_TrackletsToDIASource'
    n = dbc.execute(sql)
    n = dbc.fetchone()[0]
    if(not n):
        raise(DBSetupError('Please populate the mops_TrackletsToDIASource table.'))
    return



if(__name__ == '__main__'):
    sanityCheck()        
    unittest.main()





















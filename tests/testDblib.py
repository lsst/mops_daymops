#!/usr/bin/env python
import math
import random
import unittest

import MySQLdb as DBI

try:
    import lsst.mops.daymops.dblib as dblib
    from lsst.mops.daymops.DiaSource import DiaSource
except:
    raise(ImportError('Please setup daymops first.'))




# Constants.
DB_HOST = 'localhost'
DB_PORT = '3306'
DB_DBNAME = 'mops_test'
DB_USER = 'www'
DB_PASS = 'zxcvbnm'




class TestDblib(unittest.TestCase):
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
        
        # Retrieve the list of DiaSource instances in the database.
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
from DIASource order by DIASource.diaSourceId'''
        
        # Send the query.
        n = self._dbc.execute(sql)
        
        # Create the DiaSource objects manually.
        self.trueDiaSources = {}
        self.trueDiaSourceIds = []
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
            self.trueDiaSourceIds.append(_id)   # Just to check ordering.
            
            row = self._dbc.fetchone()
        return
    
    def testSimpleObjectFetch(self):
        # Retrieve all known DiaSources by using dblib.simpleObjectFetch()
        cols = [('diaSourceId', 'Long'),
                ('ra',          'Double'),
                ('decl',        'Double'),
                ('filterId',    'Int'),
                ('taiMidPoint', 'Double'),
                ('obsCode',     'String'),
                ('apFlux',      'Double'),
                ('apFluxErr',   'Double'),
                ('refMag',      'Double')]
        iter = dblib.simpleObjectFetch(self.dbLocStr,
                                       table='DIASource',
                                       className='DiaSource',
                                       columns=cols,
                                       where=None,
                                       orderBy='diaSourceId')
        diaSources = {}
        diaSourceIds = []
        for d in iter:
            _id = d.getDiaSourceId()
            diaSources[_id] = d
            diaSourceIds.append(_id)
        
        # Same number?
        self.failUnlessEqual(len(self.trueDiaSources), len(diaSources),
                             'the number of DiaSources is different: %d /= %d' %(len(self.trueDiaSources), len(diaSources)))
        
        # And now, one by one.
        for _id in self.trueDiaSources.keys():
            self.failUnless(_id in diaSources.keys())
            
            for attr in ('diaSourceId', 'ra', 'dec', 'filterId', 'taiMidPoint',
                         'obsCode', 'apFlux', 'apFluxErr', 'refMag'):
                self.checkObjs(self.trueDiaSources[_id], diaSources[_id], attr)
        
        # Check ordering
        for i in range(len(self.trueDiaSourceIds)):
            self.failUnlessEqual(self.trueDiaSourceIds[i], diaSourceIds[i],
                                 'Wrong ordering: ID %d /= ID %d' %(self.trueDiaSourceIds[i], diaSourceIds[i]))
        return
    
    def testSimpleObjectFetchWhere(self):
        # Retrieve all known DiaSources by using dblib.simpleObjectFetch()
        cols = [('diaSourceId', 'Long'),
                ('ra',          'Double'),
                ('decl',        'Double'),
                ('filterId',    'Int'),
                ('taiMidPoint', 'Double'),
                ('obsCode',     'String'),
                ('apFlux',      'Double'),
                ('apFluxErr',   'Double'),
                ('refMag',      'Double')]
        iter = dblib.simpleObjectFetch(self.dbLocStr,
                                       table='DIASource',
                                       className='DiaSource',
                                       columns=cols,
                                       where='diaSourceId < 1000 and ra > 10.',
                                       orderBy=None)
        diaSources = {}
        diaSourceIds = []
        for d in iter:
            _id = d.getDiaSourceId()
            diaSources[_id] = d
            diaSourceIds.append(_id)
        
        # First filter the true dict. We do not care about ordering since it is 
        # tested below.
        trueDiaSources = dict([(_id, d) for (_id, d) in self.trueDiaSources.items() \
                               if d.getDiaSourceId() < 1000 and \
                                  d.getRa() > 10.])
        
        # Same number?
        self.failUnlessEqual(len(trueDiaSources), len(diaSources),
                             'the number of DiaSources is different: %d /= %d' %(len(trueDiaSources), len(diaSources)))
        
        # And now, one by one.
        for _id in trueDiaSources.keys():
            self.failUnless(_id in diaSources.keys())
            
            for attr in ('diaSourceId', 'ra', 'dec', 'filterId', 'taiMidPoint',
                         'obsCode', 'apFlux', 'apFluxErr', 'refMag'):
                self.checkObjs(trueDiaSources[_id], diaSources[_id], attr)
        return
    
    def testSimpleObjectFetchWhereOrdered(self):
        # Retrieve all known DiaSources by using dblib.simpleObjectFetch()
        cols = [('diaSourceId', 'Long'),
                ('ra',          'Double'),
                ('decl',        'Double'),
                ('filterId',    'Int'),
                ('taiMidPoint', 'Double'),
                ('obsCode',     'String'),
                ('apFlux',      'Double'),
                ('apFluxErr',   'Double'),
                ('refMag',      'Double')]
        iter = dblib.simpleObjectFetch(self.dbLocStr,
                                       table='DIASource',
                                       className='DiaSource',
                                       columns=cols,
                                       where='diaSourceId < 1000 and ra > 10.',
                                       orderBy=None)
        diaSources = {}
        diaSourceIds = []
        for d in iter:
            _id = d.getDiaSourceId()
            diaSources[_id] = d
            diaSourceIds.append(_id)
        
        # First filter the true dict.
        trueDiaSources = {}
        trueDiaSourceIds = []
        # Keep them ordered.
        for _id in self.trueDiaSourceIds:
            d = self.trueDiaSources[_id]
            if(d.getDiaSourceId() < 1000 and d.getRa() > 10.):
                trueDiaSourceIds.append(_id)
                trueDiaSources[_id] = d
        
        # Same number?
        self.failUnlessEqual(len(trueDiaSources), len(diaSources),
                             'the number of DiaSources is different: %d /= %d' %(len(trueDiaSources), len(diaSources)))
        
        # And now, one by one.
        for _id in trueDiaSources.keys():
            self.failUnless(_id in diaSources.keys())
            
            for attr in ('diaSourceId', 'ra', 'dec', 'filterId', 'taiMidPoint',
                         'obsCode', 'apFlux', 'apFluxErr', 'refMag'):
                self.checkObjs(trueDiaSources[_id], diaSources[_id], attr)
        
        # Check ordering
        for i in range(len(trueDiaSourceIds)):
            self.failUnlessEqual(trueDiaSourceIds[i], diaSourceIds[i],
                                 'Wrong ordering: ID %d /= ID %d' %(trueDiaSourceIds[i], diaSourceIds[i]))
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





















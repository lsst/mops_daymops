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
    import lsst.mops.daymops.DiaSourceList as DiaSourceList
    from lsst.mops.daymops.DiaSource import DiaSource
except:
    raise(ImportError('Please setup daymops first.'))




# Constants.
DB_HOST = 'localhost'
DB_PORT = '3306'
DB_DBNAME = 'mops_onelunation'
DB_USER = 'jmyers'
DB_PASS = 'jmyers'




class TestDiaSourceList(unittest.TestCase):
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
        
        # Now get the IDs of the ones for tonight.
        sql = 'select DIASourceIDTonight.DIASourceId from DIASourceIDTonight'
        n = self._dbc.execute(sql)
        self.trueIdsForTonight = [int(r[0]) for r in self._dbc.fetchall()]
        
        
        # Now fetch the sources for tonight using the function we are testing.
        # Full fetch.
        iter = DiaSourceList.diaSourceListForTonight(self.dbLocStr,
                                                     sliceId=0,
                                                     numSlices=1)
        self.tonightDiaSources = dict([(d.getDiaSourceId(), d) for d in iter])
        
        # Two partial fetches.
        iter = DiaSourceList.diaSourceListForTonight(self.dbLocStr,
                                                     sliceId=0,
                                                     numSlices=2)
        self.tonightDiaSourcesParallel = \
            dict([(d.getDiaSourceId(), d) for d in iter])
        iter = DiaSourceList.diaSourceListForTonight(self.dbLocStr,
                                                     sliceId=1,
                                                     numSlices=2)
        self.tonightDiaSourcesParallel.update(dict([(d.getDiaSourceId(), d) \
                                                    for d in iter]))
        return
    
    def testDiaSourceListForTonightQuick(self):
        return(self._testDiaSourceListForTonightQuick(self.tonightDiaSources))
    
    def testDiaSourceListForTonightQuickParallel(self):
        return(self._testDiaSourceListForTonightQuick(self.tonightDiaSourcesParallel))
    
    def testDiaSourceListForTonightFull(self):
        return(self._testDiaSourceListForTonightFull(self.tonightDiaSources))
    
    def testDiaSourceListForTonightFullParallel(self):
        return(self._testDiaSourceListForTonightFull(self.tonightDiaSourcesParallel))
    
    def testComputeVelocityStatsZero(self):
        # We do not test the non zero case since it is a trivial application of
        # lib.sphericalDistance() which is tested in testLib.py.
        # Get a random id.
        max = len(self.trueDiaSources.keys()) - 1
        idx = int(max * random.random())
        _id = self.trueDiaSources.keys()[idx]
        
        d = self.trueDiaSources[_id]
        self.failUnlessRaises(Exception, 
                              DiaSourceList.computeVelocityStats, 
                              [d, d])
        return
        
    def _testDiaSourceListForTonightQuick(self, testDataDict):
        # Check that their number is correct.
        self.failUnless(len(testDataDict.keys()) == \
                        len(self.trueIdsForTonight),
                        'incorrect number of DiaSources for tonight.')
        
    def _testDiaSourceListForTonightFull(self, testDataDict):
        # Now check that their attributes match.
        for _id in self.trueIdsForTonight:
            trueSource = self.trueDiaSources.get(_id, None)
            source = testDataDict.get(_id, None)
            
            # Check for existance.
            self.failUnless(trueSource != None, 'database inconsistency error.')
            self.failUnless(source != None, 'no DiaSource for ID %d.' %(_id))
            
            # Check for equality.
            for attr in ('diaSourceId', 'ra', 'dec', 'filterId', 'taiMidPoint', 
                         'obsCode', 'apFlux', 'apFluxErr', 'refMag'):
                self.checkObjs(trueSource, source, attr)
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





















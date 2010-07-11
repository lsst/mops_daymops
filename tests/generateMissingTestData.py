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

import random
import MySQLdb as DBI


dbh = DBI.connect(db='mops_test', user='fpierfed', passwd='f1gata')
dbc = dbh.cursor()




if(False):
    sql = 'select movingObjectId from MovingObject'
    n = dbc.execute(sql)
    rows = dbc.fetchall()
    ids = zip(*rows)[0]
    n = len(ids)
    
    # We want 1000 moving objects!
    chosenIds = []
    for i in range(1000):
        _id = str(int(random.random() * n))
        while(_id in chosenIds):
            _id = str(int(random.random() * n))
        chosenIds.append(_id)
    
    # remove everything else.
    sql = 'delete from MovingObject where movingObjectId not in (%s)'
    n = dbc.execute(sql %(','.join(chosenIds)))
if(False):
    # generate fake MovingObject -> Tracklet info.
    sql = 'select movingObjectId from MovingObject'
    n = dbc.execute(sql)
    rows = dbc.fetchall()
    moIds = zip(*rows)[0]
    noTot = len(moIds)
    
    sql = 'select trackletId from mops_Tracklet'
    n = dbc.execute(sql)
    rows = dbc.fetchall()
    tIds = zip(*rows)[0]
    tTot = len(tIds)
    
    # we want each MovingObject to have between 0 and 10 tracklets. Overlaps are
    # fine.
    sql = 'insert into mops_MovingObjectToTracklet (movingObjectId, trackletId)'
    sql += ' values '
    for moId in moIds:
        valSQL = ''
        n = int(random.uniform(0, 10))
        if(not n):
            continue
        # choose n random Tracklets.
        idxs = []
        for i in range(n):
            idx = int(random.random() * (tTot - 1))
            while(idx in idxs):
                idx = int(random.random() * (tTot - 1))
            idxs.append(idx)
            valSQL += '(%d, %d), ' %(moId, tIds[idx])
        valSQL = valSQL[:-2]
        
        # Insert the values.
        # print(sql + valSQL)
        dbc.execute(sql + valSQL)
if(False):
    # Generate random Track linkages.
    sql = 'select trackletId from mops_Tracklet'
    n = dbc.execute(sql)
    rows = dbc.fetchall()
    trIds = zip(*rows)[0]
    trTot = len(trIds)
    
    # Generate 1000 trackIds
    tIds = range(1, 1001, 1)
    tTot = len(tIds)
    
    # we want each Track to have between 0 and 10 tracklets. Overlaps are fine.
    sql = 'insert into mops_TracksToTracklet (trackId, trackletId) values '
    for tId in tIds:
        valSQL = ''
        n = int(random.uniform(0, 10))
        if(not n):
            continue
        # choose n random Tracklets.
        idxs = []
        for i in range(n):
            idx = int(random.uniform(1, trTot))
            while(idx in idxs):
                idx = int(random.uniform(1, trTot))
            idxs.append(idx)
            valSQL += '(%d, %d), ' %(tId, trIds[idx])
        valSQL = valSQL[:-2]
        
        # Insert the values.
        # print(sql + valSQL)
        dbc.execute(sql + valSQL)























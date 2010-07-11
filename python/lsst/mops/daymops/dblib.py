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

"""
Generic database-related functions and classes.

Here you will find funtions that retrieve an object from the database and build
an instance of a given class based on the data retrieved.
"""
from SafeDbStorage import SafeDbStorage
import lsst.daf.persistence as persistence

# Supported classes
from MovingObject import MovingObject
from Orbit import Orbit
from Tracklet import Tracklet
from DiaSource import DiaSource

# import time





def simpleObjectFetch(dbLocStr, table, className, columns, where=None, 
                      orderBy=None):
    """
    Fetch relevant rows from a given table and instantiate one object per row.
    It is assumed that the objects to be created have a method called
        setCol_i where col_i is columns[i][0]
    Also, it is essential that one can instantiate the class without passing any
    argument to the constructor. For instance:
        obj = className()
    
    The SQL used is
        select <col1>, <col2>, <col3>[, ...] from table where <where>
    
    Column types are specified together with their names. For instance
        [('trackletId', 'Long'), ('valTot', 'Double'), ...]
    Supported types are those defined in lsst.daf.persistence.
    
    @param dbLocStr: database connection string.
    @param table: the name of the database table to select from.
    @param className: the name of the class to instantiate.
    @param columns: the list of column names, types: [(col1, type1), ...]
    @param where: the where SQL statement.
    
    Return
        Iterator: [obj1, obj2, ...]
    """
    # Simple sanity check.
    if([c for c in columns if len(c) != 2]):
        raise(Exception('columns must specify both column name and type.'))
    
    if(className not in globals().keys()):
        msg = '%s is not supported yet for simple DB extraction' % (className)
        raise(NotImplementedError(msg))
    
    # Send the query.
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    db.setTableForQuery(table)
    [db.outColumn(c[0]) for c in columns]
    if(where):
        db.setQueryWhere(where)
    if(orderBy):
        db.orderBy(orderBy)
    db.query()
    
    # Fetch the results and instantiate the objects.
    while(db.next()):
        yield(_simpleObjectCreation(db, className, columns))
    db.finishQuery()
    # return


def _simpleObjectCreation(db, name, cols):
    """
    This is a bit of black magic, sorry.
    """
    obj = globals()[name]()
    # This would be significantly faster: approximately 12% faster but we canot
    # do it since we need custom setters (especialy for classes coming from C++
    # setters = [lambda x: setattr(obj, '_%s' % (c[0]), x) for c in cols]
    setters = [getattr(obj, 'set%s%s' % (c[0][0].upper(), c[0][1:])) \
               for c in cols]
    fetchers = [getattr(db, 'getColumnByPos%s' % (c[1])) for c in cols]
    idxs = range(len(cols))
    _setAttrs(setters, fetchers, idxs)
    return(obj)


def _setAttrs(setters, fetchers, indeces):
    """
    More black magic: this basically does setYYY(getColumnByPosXXX(j)) for the
    appropriate values of j.
    """
    [setters[i](fetchers[i](indeces[i])) for i in range(len(indeces))]
    return


def simpleTwoObjectFetch(dbLocStr, table, className1, columns1, 
                         className2, columns2, where=None, orderBy=None):
    """
    Fetch relevant rows from one table and instantiate two objects per row. The
    main use case for this is
        Generate two objects/row from a single table (e.g. MovingObject and
        Orbit from the MovingObject table).
    
    It is assumed that the objects to be created have a method called
        setCol_i where col_i is columns_j[i][0] j=1, 2
    Also, it is essential that one can instantiate the class without passing any
    argument to the constructor. For instance:
        obj = className()
        
    The SQL used is
        select <col1>, <col2>, <col3>[, ...] from table where <where>
    
    Column types are specified together with their names. For instance
        [('trackletId', 'Long'), ('valTot', 'Double'), ...]
    Supported types are those defined in lsst.daf.persistence.
    
    @param dbLocStr: database connection string.
    @param table: the name of the database table to select from.
    @param className1: the name of the first class to instantiate.
    @param className2: the name of the second class to instantiate.
    @param columns1: the list of column names, types for obj1: [(col1, type1), ]
    @param columns2: the list of column names, types for obj2: [(col1, type1), ]
    @param where: the where SQL statement.
    
    Return
        Iterator: [(obj11, obj21), (obj12, obj22), ...]
    """
    # Simple sanity check.
    if([c for c in columns1+columns2 if len(c) != 2]):
        raise(Exception('columns must specify both column name and type.'))
    
    if(className1 not in globals().keys()):
        msg = '%s is not supported yet for simple DB extraction' % (className1)
        raise(NotImplementedError(msg))
    if(className2 not in globals().keys()):
        msg = '%s is not supported yet for simple DB extraction' % (className2)
        raise(NotImplementedError(msg))
    
    # Send the query.
    # t0 = time.time()
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    db.setTableForQuery(table)
    [db.outColumn(c[0]) for c in columns1+columns2]
    if(where):
        db.setQueryWhere(where)
    if(orderBy):
        db.orderBy(orderBy)
    db.query()
    # print('%.02fs compose and send query' % (time.time() - t0))
    
    # Fetch the results and instantiate the objects.
    # tt0 = time.time()
    class1 = globals()[className1]
    class2 = globals()[className2]
    setterNames1 = ['set%s%s' % (c[0][0].upper(), c[0][1:]) for c in columns1]
    setterNames2 = ['set%s%s' % (c[0][0].upper(), c[0][1:]) for c in columns2]
    fetchers1 = [getattr(db, 'getColumnByPos%s' % (c[1])) for c in columns1]
    fetchers2 = [getattr(db, 'getColumnByPos%s' % (c[1])) for c in columns2]
    
    # Compute the column indeces.
    idxs1 = range(len(columns1))
    idxs2 = range(len(columns1), len(columns1 + columns2), 1)
    
    while(db.next()):
        # t0 = time.time()
        o1 = class1()
        o2 = class2()
        setters1 = [getattr(o1, name) for name in setterNames1]
        setters2 = [getattr(o2, name) for name in setterNames2]
        
        _setAttrs(setters1, fetchers1, idxs1)
        _setAttrs(setters2, fetchers2, idxs2)
        # print('%.02fs fetch one row' % (time.time() - t0))
        yield((o1, o2))
    # print('%.02fs fetch all rows' % (time.time() - tt0))
    db.finishQuery()
    # return



def profileThis():
    """
    Utility funtion for profiling purposes only.
    """
    for moOrb in simpleTwoObjectFetch('mysql://localhost:3306/mops_onelunation',
                                      'MovingObject',
                                      'MovingObject',
                                      [('movingObjectId', 'Long'), 
                                       ('mopsStatus', 'String'),
                                       ('h_v', 'Double'), 
                                       ('g', 'Double'), 
                                       ('arcLengthDays', 'Double')],
                                      'Orbit',
                                      [('q', 'Double'), 
                                       ('e', 'Double'), 
                                       ('i', 'Double'), 
                                       ('node', 'Double'), 
                                       ('argPeri', 'Double'), 
                                       ('timePeri', 'Double'), 
                                       ('epoch', 'Double'), 
                                       ('src01', 'Double'), 
                                       ('src02', 'Double'), 
                                       ('src03', 'Double'), 
                                       ('src04', 'Double'), 
                                       ('src05', 'Double'), 
                                       ('src06', 'Double'), 
                                       ('src07', 'Double'), 
                                       ('src08', 'Double'), 
                                       ('src09', 'Double'), 
                                       ('src10', 'Double'), 
                                       ('src11', 'Double'), 
                                       ('src12', 'Double'), 
                                       ('src13', 'Double'), 
                                       ('src14', 'Double'), 
                                       ('src15', 'Double'), 
                                       ('src16', 'Double'), 
                                       ('src17', 'Double'), 
                                       ('src18', 'Double'), 
                                       ('src19', 'Double'), 
                                       ('src20', 'Double'), 
                                       ('src21', 'Double')],
                                      where='mopsStatus != "M" and ' + \
                                            'movingObjectId < 50000'):
        # print(moOrb)
        pass
    return



if(__name__ == '__main__'):
    import cProfile
    import pstats
    
    
    f = '/tmp/daymops_profile'
    cProfile.run('profileThis()', f)
    p = pstats.Stats(f)
    p.strip_dirs()
    p.sort_stats('cumulative')
    p.print_stats()











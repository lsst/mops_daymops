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
Helper functions to create/retrieve a list of Track instances.
"""
from DayMOPSObject import DayMOPSObject
from DiaSource import DiaSource
from Tracklet import Tracklet
import DiaSourceList
from Track import Track
from SafeDbStorage import SafeDbStorage

import lsst.daf.persistence as persistence




def newTracks(dbLocStr, shallow=True, sliceId=None, numSlices=None):
    """
    Fetch all Tracks.
    
    Use  sliceId and numSlices to implement some form of parallelism.
    If shallow=False, then fetch the DIASources also.
    
    Remember that right now we are using the mops_TracksToTracklet as a 
    temporary table as there is no need to store Tracks permanently. This is why
    we can get away with this simple strategy.
    
    @param dbLocStr: database connection string.
    @param shallow: if True, do not bother retrieving DIASources per Tracklet.
    @param sliceId: Id of the current Slice.
    @param numSlices: number of available slices (i.e. MPI universe size - 1)
    
    Return
    Interator to the list of Track instances.
    """
    return(_fetchTracks(dbLocStr, '', [], shallow, sliceId, numSlices))


def _fetchTracks(dbLocStr, where, extraTables=[], shallow=True, 
                 sliceId=None, numSlices=None):
    """
    Fetch all Tracks.
    
    The teables from which we select (in ddition to mops_TracksToTracklet) have 
    to be specified in the extraTables list.
    
    The SQL where clause has to be specified. In it, specify the full 
    table names (e.g. mops_Track.status instead of just status).
    
    Use  sliceId and numSlices to implement some form of parallelism.
    If shallow=False, then fetch the DIASources also.
    
    Return interator.
    """
    if(shallow):
        return(_fetchShallowTracks(dbLocStr, where, extraTables, sliceId, 
                                   numSlices))
    # else:
    return(_fetchDeepTracks(dbLocStr, where, extraTables, sliceId, 
                            numSlices))




def _fetchShallowTracks(dbLocStr, where, extraTables=[], sliceId=None, 
                        numSlices=None):
    """
    Fetch all Tracks.
    
    The teables from which we select (in ddition to mops_TracksToTracklet) have 
    to be specified in the extraTables list.
    
    The SQL where clause has to be specified. In it, specify the full 
    table names (e.g. mops_Track.status instead of just status).
    
    Use  sliceId and numSlices to implement some form of parallelism.
    
    Do not fetch DiaSources (shallow fetch).
    
    Return interator.
    """
    # Send the query.
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    tables = ['mops_TracksToTracklet', ] + list(extraTables)
    db.setTableListForQuery(tables)
    db.outColumn('mops_TracksToTracklet.trackId')
    db.outColumn('mops_TracksToTracklet.trackletId')
    
    w = ''
    if(where):
        w = where
    if(sliceId != None and numSlices > 1):
        w += ' and mops_TracksToTracklet.trackId %% %d = %d' \
             %(numSlices, sliceId)
    db.setQueryWhere(w)
    db.orderBy('mops_TracksToTracklet.trackId')
    db.query()
    
    # Fetch the results.
    refId = None
    t = Track()
    tracklets = []
    while(db.next()):
        # DiaSource
        tr = Tracklet()
        tr.setTrackletId(db.getColumnByPosLong(1))
        
        trackId = db.getColumnByPosLong(0)
        if(refId == None):
            # First pass.
            t.setTrackId(trackId)
            
            # Init refId.
            refId = trackId
        elif(refId != trackId):
            # New Track. Finish up the old track.
            t.setTracklets(tracklets)
            yield(t)
            
            # Now reset the diaSources list.
            tracklets = []
            
            # And reset t.
            t = Track()
            t.setTrackId(trackId)
            
            # Update refId.
            refId = trackId
        else:
            # Same track: do not modify it.
            pass
        
        # Add tr to tracklets.
        tracklets.append(tr)
    db.finishQuery()
    
    # yield the last one.
    if(tracklets):
        t.setTracklets(tracklets)
        yield(t)
    
    del(db)
    # return


def _fetchDeepTracks(dbLocStr, where, extraTables=[], sliceId=None, 
                     numSlices=None):
    """
    Fetch all Tracks.
    
    The teables from which we select (in ddition to mops_TracksToTracklet) have 
    to be specified in the extraTables list.
    
    The SQL where clause has to be specified. In it, specify the full 
    table names (e.g. mops_Track.status instead of just status).
    
    Use  sliceId and numSlices to implement some form of parallelism.
    
    Fetch DiaSources per Tracklet (deep fetch).
    
    Return interator.
    """
    # Send the query.
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    tables = ['mops_TracksToTracklet', 
              'mops_TrackletsToDIASource', 
              'mops_Tracklet',
              'DiaSource'] + list(extraTables)
    db.setTableListForQuery(tables)
    db.outColumn('mops_TracksToTracklet.trackId')                       # 0
    db.outColumn('mops_Tracklet.trackletId')                            # 1
    db.outColumn('mops_Tracklet.velRa')                                 # 2
    db.outColumn('mops_Tracklet.velDecl')                               # 3
    db.outColumn('mops_Tracklet.velTot')                                # 4
    db.outColumn('mops_Tracklet.status')                                # 5
    db.outColumn('DIASource.diaSourceId')                               # 6
    db.outColumn('DIASource.ra')                                        # 7
    db.outColumn('DIASource.decl')                                      # 8
    db.outColumn('DIASource.filterId')                                  # 9
    db.outColumn('DIASource.taiMidPoint')                               # 10
    db.outColumn('DIASource.obsCode')                                   # 11
    db.outColumn('DIASource.apFlux')                                    # 12
    db.outColumn('DIASource.apFluxErr')                                 # 13
    db.outColumn('DIASource.refMag')                                    # 14
    
    w2 = 'mops_TracksToTracklet.trackletId=mops_Tracklet.trackletId'
    w2 += ' and mops_Tracklet.trackletId=mops_TrackletsToDIASource.trackletId'
    w2 += ' and mops_TrackletsToDIASource.diaSourceId=DIASource.diaSourceId'
    
    w = ''
    if(where):
        w = '%s and %s' %(w2, where)
    else:
        w = w2
    if(sliceId != None and numSlices > 1):
        w += ' and mops_TracksToTracklet.trackId %% %d = %d' %(numSlices, sliceId)
    db.setQueryWhere(w)
    db.orderBy('mops_TracksToTracklet.trackId,mops_TracksToTracklet.trackletId,DIASource.diaSourceId')
    db.query()
    
    # Fetch the results.
    refTrackId = None
    refTrackletId = None
    # create an empty Track and tracklet: we will update them as we go.
    t = Track()
    tr = Tracklet()
    tracklets = []
    diaSources = []
    while(db.next()):
        # Get the current trackId and trackletId.
        trackId = db.getColumnByPosLong(0)
        trackletId = db.getColumnByPosLong(1)
        
        # Create a new DiaSource, since we have a new one per row.
        d = DiaSource()
        d.setDiaSourceId(db.getColumnByPosLong(6))
        d.setRa(db.getColumnByPosDouble(7))
        d.setDec(db.getColumnByPosDouble(8))
        d.setFilterId(db.getColumnByPosInt(9))
        d.setTaiMidPoint(db.getColumnByPosDouble(10))
        d.setObsCode(db.getColumnByPosString(11))
        d.setApFlux(db.getColumnByPosDouble(12))
        d.setApFluxErr(db.getColumnByPosDouble(13))
        d.setRefMag(db.getColumnByPosDouble(14))
        
        # Now see if we are within the same Track and/or same Tracklet.
        if(refTrackId == None):
            # First pass: nothing is set at this point.
            t.setTrackId(trackId)
            
            tr.setTrackletId(trackletId)
            tr.setVelRa(db.getColumnByPosDouble(2))
            tr.setVelDec(db.getColumnByPosDouble(3))
            tr.setVelTot(db.getColumnByPosDouble(4))
            tr.setStatus(db.getColumnByPosString(5))
            
            # Init refTrackId and refTrackletId.
            refTrackId = trackId
            refTrackletId = trackletId
        elif(refTrackId != trackId):
            # New Track. Finish up the old track. This also means that the
            # current Tracklet needs to be closed and then re-initialized.
            tr.setDiaSources(diaSources)
            t.setTracklets(tracklets)
            yield(t)
            
            # Now reset the lists.
            diaSources = []
            tracklets = []
            
            # And reset t and tr.
            t = Track()
            t.setTrackId(trackId)
            
            tr = Tracklet()
            tr.setTrackletId(trackletId)
            tr.setVelRa(db.getColumnByPosDouble(2))
            tr.setVelDec(db.getColumnByPosDouble(3))
            tr.setVelTot(db.getColumnByPosDouble(4))
            tr.setStatus(db.getColumnByPosString(5))
            
            # Update refTrackId.
            refTrackId = trackId
            refTrackletId = trackletId
        else:
            # Same track: do not modify it. But let's see what is up with the 
            # Tracklet. We only have two options: either we are within the same 
            # Tracklet (in which case we do not do anything) or we are in a new 
            # Tracklet (in which case we update its diaSources).
            if(refTrackletId != trackletId):
                # New Tracklet!
                tr.setDiaSources(diaSources)
                
                # Now reset the diaSources list.
                diaSources = []
                
                # And reset tr.
                tr = Tracklet()
                tr.setTrackletId(trackletId)
                tr.setVelRa(db.getColumnByPosDouble(2))
                tr.setVelDec(db.getColumnByPosDouble(3))
                tr.setVelTot(db.getColumnByPosDouble(4))
                tr.setStatus(db.getColumnByPosString(5))
                
                # Update refTrackId.
                refTrackletId = trackletId
            else:
                # Same Tracklet: do not bother adding the tracklet to tracklets.
                # Just update diaSources and continue.
                diaSources.append(d)
                continue
        
        # Update the lists.
        diaSources.append(d)
        tracklets.append(tr)
    db.finishQuery()
    
    # yield the last one.
    if(tracklets):
        tracklets[-1].setDiaSources(diaSources)
        t.setTracklets(tracklets)
        yield(t)
    del(db)
    # return


def deleteAllTracks(dbLocStr):
    """
    Since Tracks are ephemeral, delete all Tracks urrently stored. We treat the
    Track table as temporary storage only.
    
    @param dbLocStr: database connection string.
    
    Return
    None
    """
    # Connect to the database.
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    db.startTransaction()
    db.executeSql('delete from mops_TracksToTracklet')
    db.endTransaction()
    return


def save(dbLocStr, tracks):
    """
    Save the input list of Track instances to the database.
    
    @param dbLocStr: database connection string.
    @param tracks: a list of Track instances.
    
    Return
    None
    """
    # Get the next available trackId.
    newTrackId = _getNextTrackId(dbLocStr)

    # Connect to the database.
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    
    # Prepare for insert.
    db.setTableForInsert('mops_TracksToTracklet')
    
    db.startTransaction()
    for track in tracks:
        # If the track has an id already, use that, otherwise use a new 
        # one.
        trackId = track.getTrackId()
        if(trackId == None):
            trackId = newTrackId
            track.setTrackId(newTrackId)
            newTrackId += 1
        
        # Insert the trackId <-> trackletIds info.
        for tracklet in track.getTracklets():
            db.setColumnLong('trackId', trackId)
            db.setColumnLong('trackletId', tracklet.getTrackletId())
            db.insertRow()
    db.endTransaction()
    return


def _getNextTrackId(dbLocStr):
    # Since trackId could be turned into an autoincrement column, it should not 
    # be 0. It will get incremented by 1 at the end of the function call...
    trackId = 0
        
    # Connect to the database.
    db = SafeDbStorage()
    db.setRetrieveLocation(persistence.LogicalLocation(dbLocStr))
    
    db.setTableForQuery('mops_TracksToTracklet')
    db.outColumn('max(trackId)', True)      # isExpr=True
    
    db.query()
    if(db.next() and not db.columnIsNull(0)):
        trackId = db.getColumnByPosLong(0)
    db.finishQuery()
    return(trackId + 1)





"""
Helper functions to create/retrieve a list of Tracklet instances.
"""
from DayMOPSObject import DayMOPSObject
from DiaSource import DiaSource
import DiaSourceList
from Tracklet import Tracklet, STATUS
import lsst.daf.persistence as persistence
from SafeDbStorage import SafeDbStorage



def newTrackletsFromTonight(dbLocStr, shallow=True, 
                            sliceId=None, numSlices=None):
    """
    Fetch unattributed Tracklets with at least one DiaSource form tonight.
    
    Use  sliceId and numSlices to implement some form of parallelism.
    If shallow=False, then fetch the DIASources also.
    
    @param dbLocStr: database connection string.
    @param shallow: if True, do not bother retrieving DIASources per Tracklet.
    @param sliceId: Id of the current Slice.
    @param numSlices: number of available slices (i.e. MPI universe size - 1)
    
    Return 
    Interator to the list of Tracklet instances.
    """
    where='mops_TrackletsToDIASource.diaSourceId=DIASourceIDTonight.DIASourceId'
    return(_fetchTracklets(dbLocStr, 
                           where,
                           ('DIASourceIDTonight', ),
                           shallow, 
                           sliceId, 
                           numSlices))


def newTracklets(dbLocStr, fromMjd=None, toMjd=None, shallow=True, 
                 sliceId=None, numSlices=None):
    """
    Fetch non attrinuted/linked/precovered (i.e. status='U') Tracklets.
    If fromMjd != None, then only consider tracklets which have at least one 
    DiaSource with taiMidPoint >= fromMjd.
    
    If toMjd != None, then only consider tracklets which have at least one 
    DiaSource with taiMidPoint <= toMjd.
    
    Use  sliceId and numSlices to implement some form of parallelism.
    If shallow=False, then fetch the DIASources also.
    
    @param dbLocStr: database connection string.
    @param shallow: if True, do not bother retrieving DIASources per Tracklet.
    @param fromMjd: min DIASource MJD to consider for Tracklet retrieval.
    @param toMjd: max DIASource MJD to consider for Tracklet retrieval.
    @param sliceId: Id of the current Slice.
    @param numSlices: number of available slices (i.e. MPI universe size - 1)
    
    Return
    Interator to the list of Tracklet instances.
    """
    where = ''
    if(fromMjd != None):
        where += 'DIASource.taiMidPoint >= %f' %(fromMjd)
    if(toMjd != None):
        if(where):
            where += 'and DIASource.taiMidPoint <= %f' %(toMjd)
        else:
            where = 'DIASource.taiMidPoint <= %f' %(toMjd)
    return(_fetchTracklets(dbLocStr, 
                           where,
                           [],
                           shallow, 
                           sliceId, 
                           numSlices))


def allTrackletsForMovingObject(dbLocStr, movingObjectId, shallow=True, 
                                sliceId=None, numSlices=None):
    """
    Fetch all Tracklet instances associated to a given MovingObject (via its
    movingObjectId).
    
    Use  sliceId and numSlices to implement some form of parallelism.
    If shallow=False, then fetch the DIASources also.
    
    @param dbLocStr: database connection string.
    @param shallow: if True, do not bother retrieving DIASources per Tracklet.
    @param movingObjectId: ID of the MovingObject associated to the Tracklets.
    @param sliceId: Id of the current Slice.
    @param numSlices: number of available slices (i.e. MPI universe size - 1)
    
    Return
    Interator to the list of Tracklet instances.
    """
    w = 'mops_MovingObjectToTracklet.movingObjectId=%d ' %(movingObjectId)
    w += 'and mops_MovingObjectToTracklet.trackletId=mops_Tracklet.trackletId'
    extra = ['mops_MovingObjectToTracklet', ]
    return(_fetchTracklets(dbLocStr, w, extra, shallow, sliceId, numSlices))


def _fetchTracklets(dbLocStr, where, extraTables=[], shallow=True, 
                    sliceId=None, numSlices=None):
    """
    Fetch non attrinuted/linked/precovered (i.e. status=STATUS['UNATTRIBUTED']) 
    Tracklets.
    
    The teables from which we select (in ddition to mops_Tracklet and 
    mops_TrackletsToDIASource) have to be specified in the extraTables list.
    
    The SQL where clause has to be specified. In it, specify the full 
    table names (e.g. mops_Tracklet.status instead of just status).
    
    Use  sliceId and numSlices to implement some form of parallelism.
    If shallow=False, then fetch the DIASources also.
    
    Return interator.
    """
    if(shallow):
        return(_fetchShallowTracklets(dbLocStr, where, extraTables, sliceId, 
                                      numSlices))
    # else:
    return(_fetchDeepTracklets(dbLocStr, where, extraTables, sliceId, 
                               numSlices))




def _fetchShallowTracklets(dbLocStr, where, extraTables=[], sliceId=None, 
                           numSlices=None):
    """
    Fetch non attrinuted/linked/precovered (i.e. status=STATUS['UNATTRIBUTED']) 
    Tracklets.
    
    The teables from which we select (in ddition to mops_Tracklet and 
    mops_TrackletsToDIASource) have to be specified in the extraTables list.
    
    The SQL where clause has to be specified. In it, specify the full 
    table names (e.g. mops_Tracklet.status instead of just status).
    
    Use  sliceId and numSlices to implement some form of parallelism.
    
    Do not fetch DiaSources (shallow fetch).
    
    Tracklets ordered by trackletId
    
    Return interator.
    """
    # Send the query.
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    tables = ['mops_Tracklet', 'mops_TrackletsToDIASource'] + list(extraTables)
    db.setTableListForQuery(tables)
    db.outColumn('distinct(mops_Tracklet.trackletId)', True)
    db.outColumn('mops_Tracklet.velRa')
    db.outColumn('mops_Tracklet.velDecl')
    db.outColumn('mops_Tracklet.velTot')
    db.outColumn('mops_Tracklet.status')
    
    w2 = 'mops_Tracklet.status="%s"' %(STATUS['UNATTRIBUTED'])
    w2 += ' and mops_Tracklet.trackletId=mops_TrackletsToDIASource.trackletId'
    
    w = ''
    if(where):
        w = '%s and %s' %(w2, where)
    else:
        w = w2
    if(sliceId != None and numSlices > 1):
        w += ' and mops_Tracklet.trackletId %% %d = %d' %(numSlices, sliceId)
    db.setQueryWhere(w)
    db.orderBy('mops_Tracklet.trackletId')
    db.query()
    
    # Fetch the results.
    while(db.next()):
        t = Tracklet()
        t.setTrackletId(db.getColumnByPosLong(0))
        t.setVelRa(db.getColumnByPosDouble(1))
        t.setVelDec(db.getColumnByPosDouble(2))
        t.setVelTot(db.getColumnByPosDouble(3))
        t.setStatus(db.getColumnByPosString(4))
        yield(t)
    db.finishQuery()
    del(db)
    # return


def _fetchDeepTracklets(dbLocStr, where, extraTables=[], sliceId=None, 
                        numSlices=None):
    """
    Fetch non attrinuted/linked/precovered (i.e. status=STATUS['UNATTRIBUTED']) 
    Tracklets.
    
    The teables from which we select (in ddition to mops_Tracklet and 
    mops_TrackletsToDIASource) have to be specified in the extraTables list.
    
    The SQL where clause has to be specified. In it, specify the full 
    table names (e.g. mops_Tracklet.status instead of just status).
    
    Use  sliceId and numSlices to implement some form of parallelism.
    
    Fetch DiaSources (deep fetch).
    
    Tracklets ordered by trackletId
    
    Return interator.
    """
    # Send the query.
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    tables = ['mops_Tracklet', 'mops_TrackletsToDIASource', 'DIASource'] + \
             list(extraTables)
    db.setTableListForQuery(tables)
    db.outColumn('mops_Tracklet.trackletId')
    db.outColumn('mops_Tracklet.velRa')
    db.outColumn('mops_Tracklet.velDecl')
    db.outColumn('mops_Tracklet.velTot')
    db.outColumn('mops_Tracklet.status')
    db.outColumn('DIASource.diaSourceId')
    db.outColumn('DIASource.ra')
    db.outColumn('DIASource.decl')
    db.outColumn('DIASource.filterId')
    db.outColumn('DIASource.taiMidPoint')
    db.outColumn('DIASource.obsCode')
    db.outColumn('DIASource.apFlux')
    db.outColumn('DIASource.apFluxErr')
    db.outColumn('DIASource.refMag')
    
    w2 = 'mops_Tracklet.status="%s"' %(STATUS['UNATTRIBUTED'])
    w2 += ' and mops_Tracklet.trackletId=mops_TrackletsToDIASource.trackletId'
    w2 += ' and DIASource.diaSourceId=mops_TrackletsToDIASource.diaSourceId'
    
    w = ''
    if(where):
        w = '%s and %s' %(w2, where)
    else:
        w = w2
    if(sliceId != None and numSlices > 1):
        w += ' and mops_Tracklet.trackletId %% %d = %d' %(numSlices, sliceId)
    db.setQueryWhere(w)
    db.orderBy('mops_Tracklet.trackletId')
    db.query()
    
    # Fetch the results.
    refId = None
    t = Tracklet()
    diaSources = []
    while(db.next()):
        # DiaSource
        d = DiaSource()
        d.setDiaSourceId(db.getColumnByPosLong(5))
        d.setRa(db.getColumnByPosDouble(6))
        d.setDec(db.getColumnByPosDouble(7))
        d.setFilterId(db.getColumnByPosInt(8))
        d.setTaiMidPoint(db.getColumnByPosDouble(9))
        d.setObsCode(db.getColumnByPosString(10))
        d.setApFlux(db.getColumnByPosDouble(11))
        d.setApFluxErr(db.getColumnByPosDouble(12))
        d.setRefMag(db.getColumnByPosDouble(13))
        
        trackletId = db.getColumnByPosLong(0)
        if(refId == None):
            # First pass.
            t.setTrackletId(trackletId)
            t.setVelRa(db.getColumnByPosDouble(1))
            t.setVelDec(db.getColumnByPosDouble(2))
            t.setVelTot(db.getColumnByPosDouble(3))
            t.setStatus(db.getColumnByPosString(4))
            
            # Init refId.
            refId = trackletId
        elif(refId != trackletId):
            # New Tracklet. Finish up theold tracklet.
            t.setDiaSources(diaSources)
            yield(t)
            
            # Now reset the diaSources list.
            diaSources = []
            
            # And reset t.
            t = Tracklet()
            t.setTrackletId(trackletId)
            t.setVelRa(db.getColumnByPosDouble(1))
            t.setVelDec(db.getColumnByPosDouble(2))
            t.setVelTot(db.getColumnByPosDouble(3))
            t.setStatus(db.getColumnByPosString(4))
            
            # Update refId.
            refId = trackletId
        else:
            # Same tracklet: do not modify it.
            pass
        
        # Add d to diaSources.
        diaSources.append(d)
    db.finishQuery()
    
    # yiled the last one.
    if(diaSources):
        t.setDiaSources(diaSources)
        yield(t)
    
    del(db)
    # return


def getArcLength(tracklets):
    """
    Compute the arc length in days from the DiaSources that are part of this
    list of Tracklets.
    
    @param tracklets: a list of Tracklet instances.
    
    Return
    Total arc length in days for all the DIASources part of the unition of input
    Tracklet instances.
    """
    # TODO: is this the most efficient way of computing the arc length?
    # Probably not, but at least is pretty elegant.
    diaSources = reduce(lambda a, b: a+b, 
                        [t.getDiaSources() for t in tracklets])
    return(DiaSourceList.getTimeSpan(diaSources))


def save(dbLocStr, tracklets):
    """
    Save the input list of Tracklet instances to the database.
    
    @param dbLocStr: database connection string.
    @param tracklets: a list of Tracklet instances.
    
    Return
    None
    """
    # Get the next available trackletId.
    newTrackletId = _getNextTrackletId(dbLocStr)

    # Connect to the database.
    dbTrk = SafeDbStorage()
    dbTrk.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    dbSrc = SafeDbStorage()
    dbSrc.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    
    # Prepare for insert.
    dbTrk.setTableForInsert('mops_Tracklet')
    dbSrc.setTableForInsert('mops_TrackletsToDIASource')
    
    # dbTrk.startTransaction()
    # dbSrc.startTransaction()
    for tracklet in tracklets:
        # If the tracklet has an id already, use that, otherwise use a new 
        # one.
        trackletId = tracklet.getTrackletId()
        velRa = tracklet.getVelRa()
        velDec = tracklet.getVelDec()
        velTot = tracklet.getVelTot()
        if(trackletId == None):
            trackletId = newTrackletId
            tracklet.setTrackletId(newTrackletId)
            newTrackletId += 1
        
        # Insert values.
        dbTrk.setColumnLong('trackletId', trackletId)
        dbTrk.setColumnString('status', tracklet.getStatus())
        if(velRa != None and velDec != None and velTot != None):
            dbTrk.setColumnDouble('velRa', velRa)
            dbTrk.setColumnDouble('velDecl', velDec)
            dbTrk.setColumnDouble('velTot', velTot)
        else:
            dbTrk.setColumnToNull('velRa')
            dbTrk.setColumnToNull('velDecl')
            dbTrk.setColumnToNull('velTot')
        dbTrk.insertRow()
        
        # Insert the trackletId <-> diaSourceIds info.
        for diaSource in tracklet.getDiaSources():
            dbSrc.setColumnLong('trackletId', trackletId)
            dbSrc.setColumnLong('diaSourceId', diaSource.getDiaSourceId())
            dbSrc.insertRow()
    # dbTrk.endTransaction()
    # dbSrc.endTransaction()
    return


def update(dbLocStr, tracklets):
    """
    Updates the input list of Tracklet instances in the database. It assumes
    that these Tracklets are already present in the DB and hence only does an
    update. The non obvious implication is that the Tracklets have to have valid
    trackletIds. It only updates the Tracklets themselves and their associations
    to DiaSources, not the DiaSource themselves. For that, update the DiaSources
    separately.
    
    @param dbLocStr: database connection string.
    @param tracklets: a list of Tracklet instances.
    
    Return
    None
    """
    # Connect to the database.
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    
    # Prepare the SQL statements.
    trSql = 'update mops_Tracklet set status="%s", velRa=%s, velDecl=%s, '
    trSql += 'velTot=%s where trackletId=%d'
    
    tdSql = 'update mops_TrackletsToDIASource set diaSourceId=%d '
    tdSql += 'where trackletId=%d'
    
    db.startTransaction()
    for tracklet in tracklets:
        # Update Tracklet data.
        trackletId = tracklet.getTrackletId()
        if(trackletId == None):
            raise(Exception('Tracklet instance %s does not have a valid ID.' \
                            %(str(tracklet))))
        velRa = tracklet.getVelRa()
        if(velRa == None):
            valRa = 'NULL'
        velDec = tracklet.getVelDec()
        if(velDec == None):
            velDec = 'NULL'
        velTot = tracklet.getVelTot()
        if(velTot == None):
            velTot = 'NULL'
        
        # Update the Tracklet data.
        db.executeSql(trSql %(tracklet.getStatus(), str(velRa), str(velDec), 
                              str(velTot), trackletId))
                
        # Update the trackletId <-> diaSourceIds info.
        for diaSource in tracklet.getDiaSources():
            db.executeSql(tdSql %(diaSource.getDiaSourceId(), trackletId))
    db.endTransaction()
    return


def updateStatus(dbLocStr, tracklets, inTransaction=False):
    """
    Updates the input list of Tracklet instances in the database. It assumes
    that these Tracklets are already present in the DB and hence only does an
    update. The non obvious implication is that the Tracklets have to have valid
    trackletIds. It only updates the Tracklets status.
    
    This is a convenience function that does a fraction of what update() does.
    However if one only needs to update status, this is faster.
    
    If inTransaction=True, do not start a new transaction since we are already 
    inside one.
    
    @param dbLocStr: database connection string.
    @param tracklets: a list of Tracklet instances.
    @param inTransaction: are we already inside a transaction?
    
    Return
    None
    """
    # Connect to the database.
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    
    # Prepare the SQL statements.
    trSql = 'update mops_Tracklet set status="%s" where trackletId=%d'
    
    if(not inTransaction):
        db.startTransaction()
    for tracklet in tracklets:
        # Update Tracklet data.
        trackletId = tracklet.getTrackletId()
        if(trackletId == None):
            raise(Exception('Tracklet instance %s does not have a valid ID.' \
                            %(str(tracklet))))
        
        # Update the Tracklet data.
        db.executeSql(trSql %(tracklet.getStatus(), trackletId))
    if(not inTransaction):
        db.endTransaction()
    return


def _getNextTrackletId(dbLocStr):
    # Since trackletId is autoincrement, it cannot be 0. It will get
    # incremented by 1 at the end of the function call...
    trackletId = 0
        
    # Connect to the database.
    db = SafeDbStorage()
    db.setRetrieveLocation(persistence.LogicalLocation(dbLocStr))
    
    db.setTableForQuery('mops_Tracklet')
    db.outColumn('max(trackletId)', True)      # isExpr=True
    
    db.query()
    if(db.next() and not db.columnIsNull(0)):
        trackletId = db.getColumnByPosLong(0)
    db.finishQuery()
    return(trackletId + 1)





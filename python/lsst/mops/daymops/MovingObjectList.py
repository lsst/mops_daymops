"""
Helper functions to create/retrieve a list of MovingObject instances.
"""
from DayMOPSObject import DayMOPSObject
from MovingObject import MovingObject, STATUS
from Orbit import STABLE_STATUS
import TrackletList
import dblib

import lsst.daf.persistence as persistence
from SafeDbStorage import SafeDbStorage




def getAllMovingObjects(dbLocStr, shallow=True, sliceId=None, numSlices=None):
    """
    Fetch all active and non merged MovingObjects we know anything about.
    
    Use  sliceId and numSlices to implement some form of parallelism.
    If shallow=False, then fetch the Tracklets also.
    
    @param dbLocStr: database connection string.
    @param shallow: if True, do not bother retrieving Tracklets per MovingObject
    @param sliceId: Id of the current Slice.
    @param numSlices: number of available slices (i.e. MPI universe size - 1)
    
    Return 
    Interator to the list of MovingObject instances.
    """
    return(_getMovingObjects(dbLocStr, 'mopsStatus != "%s"' %(STATUS['MERGED']),
                             shallow, sliceId, numSlices))

def getAllUnstableMovingObjects(dbLocStr, shallow=True, sliceId=None, 
                                numSlices=None):
    """
    Fetch all active, non merged, non stable MovingObjects we know anything 
    about.
    
    Use  sliceId and numSlices to implement some form of parallelism.
    If shallow=False, then fetch the Tracklets also.
    
    @param dbLocStr: database connection string.
    @param shallow: if True, do not bother retrieving Tracklets per MovingObject
    @param sliceId: Id of the current Slice.
    @param numSlices: number of available slices (i.e. MPI universe size - 1)
    
    Return 
    Interator to the list of MovingObject instances.
    """
    where = 'mopsStatus != "%s" and stablePass != "%s"' \
            %(STATUS['MERGED'], STABLE_STATUS['STABLE'])
    return(_getMovingObjects(dbLocStr, where, shallow, sliceId, numSlices))


def getAllPreliminaryMovingObjects(dbLocStr, shallow=True, sliceId=None, 
                                   numSlices=None):
    """
    Fetch all active, non merged, preliminary MovingObjects we know anything 
    about.
    
    Use  sliceId and numSlices to implement some form of parallelism.
    If shallow=False, then fetch the Tracklets also.
    
    @param dbLocStr: database connection string.
    @param shallow: if True, do not bother retrieving Tracklets per MovingObject
    @param sliceId: Id of the current Slice.
    @param numSlices: number of available slices (i.e. MPI universe size - 1)
    
    Return 
    Interator to the list of MovingObject instances.
    """
    where = 'mopsStatus = "%s"' %(STATUS['PRELIMINARY'])
    return(_getMovingObjects(dbLocStr, where, shallow, sliceId, numSlices))


def _getMovingObjects(dbLocStr, where, shallow=True, 
                      sliceId=None, numSlices=None):
    """
    Fetch all MovingObjects specifying the SQL where clause we want to use.
    
    Use  sliceId and numSlices to implement some form of parallelism.
    If shallow=False, then fetch the Tracklets also.
    
    Return an iterator.
    """
    if(not shallow):
        return(_getMovingObjectsDeep(dbLocStr,where,shallow,sliceId,numSlices))
    return(_getMovingObjectsShallow(dbLocStr,where,shallow,sliceId,numSlices))


def _getMovingObjectsDeep(dbLocStr, where, shallow=True, 
                          sliceId=None, numSlices=None):
    # What we do here is do a shallow fetch first and then loop over and get the
    # Tracklets.
    movingObjectsIter = _getMovingObjectsShallow(dbLocStr,
                                                 where,
                                                 shallow,
                                                 sliceId,
                                                 numSlices)
    for movingObject in movingObjectsIter:
        _id = movingObject.getMovingObjectId()
        trackletIter = TrackletList.allTrackletsForMovingObject(dbLocStr,
                                                                _id,
                                                                False,
                                                                None,
                                                                None)
        
        tracklets = [t for t in trackletIter]
        movingObject.setTracklets(tracklets)
        yield(movingObject)
    # return


def _getMovingObjectsShallow(dbLocStr, where, shallow=True, 
                             sliceId=None, numSlices=None):
    # Send the query.
    # sql: select movingObjectId, mopsStatus, h_v, g, 
    #      q, e, i, node, argPeri, timePeri, epoch, 
    #      src01, src02, src03, src04, src05, src06, src07, src08, src09, 
    #      src10, src11, src12, src13, src14, src15, src16, src17, src18, 
    #      src19, src20, src21 from MovingObject where
    #      mopsStatus != "M"
    cols = [('movingObjectId', 'Long'), 
            ('mopsStatus', 'String'),
            ('h_v', 'Double'), 
            ('g', 'Double'), 
            ('arcLengthDays', 'Double'),
            ('q', 'Double'), 
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
            ('src21', 'Double'),
            ('stablePass', 'String')]
    if(sliceId != None and numSlices > 1):
        where += ' and movingObjectId %% %d = %d' %(numSlices, sliceId)
    
    # Fetch MovingObjects and their Orbit.
    for (mo, o) in dblib.simpleTwoObjectFetch(dbLocStr,
                                              table='MovingObject',
                                              className1='MovingObject',
                                              columns1=cols[:5],
                                              className2='Orbit',
                                              columns2=cols[5:],
                                              where=where):
        # Patch the src
        o.setSrc([getattr(o, 'getSrc%02d' %(i))() for i in range(1, 22, 1)])
        # print([getattr(o, 'getSrc%02d' %(i))() for i in range(1, 22, 1)])
        
        # Now add the Orbit to each MovingObject.
        mo.setOrbit(o)
        yield(mo)
    # return


def updateStatus(dbLocStr, movingObjects, updateTrackletStatus=True):
    """
    Update the mopsStatus field of the input list of MovingObject instances in 
    the database. It is assumed that the input MovingObjects have already been
    saved in the database, of course. If not, use save().
    
    @param dbLocStr: database connection string.
    @param movingObjects: a list of MovingObject instances.
    @param updateTrackletStatus: if True, also update Tracklet status as a 
           convenience. It is assumed that the caller has properly updated the
           status of the relevant Tracklet instances.
    
    Return
    None
    """
    # Connect to the database.
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    
    # Prepare for the SQL statement.
    sql = 'update MovingObject set mopsStatus="%s" where movingObjectId=%d'
    
    db.startTransaction()
    for movingObject in movingObjects:
        movingObjectId = movingObject.getMovingObjectId()
        status = movingObject.getStatus()
        
        # Update the MovingObject data.
        db.executeSql(sql %(status, movingObjectId))
        
        # If requested, update the Tracklet status as well.
        if(updateTrackletStatus):
            TrackletList.updateStatus(dbLocStr,
                                      movingObject.getTracklets(),
                                      inTransaction=True)
    db.endTransaction()
    return



def save(dbLocStr, movingObjects, updateTrackletStatus=True):
    """
    Save the input list of MovingObject instances to the database.
    
    @param dbLocStr: database connection string.
    @param movingObjects: a list of MovingObject instances.
    @param updateTrackletStatus: if True, also update Tracklet status as a 
           convenience. It is assumed that the caller has properly updated the
           status of the relevant Tracklet instances.
    
    Return
    None
    """
    # Get the next available movingObjectId.
    newMovingObjectId = _getNextMovingObjectId(dbLocStr)

    # Connect to the database.
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    
    # Prepare for insert.
    db.setTableForInsert('MovingObject')
    
    for movingObject in movingObjects:
        # If the movingObject has an id already, use that, otherwise use a new 
        # one.
        movingObjectId = movingObject.getMovingObjectId()
        if(movingObjectId == None):
            movingObjectId = newMovingObjectId
            movingObject.setMovingObjectId(newMovingObjectId)
            newMovingObjectId += 1
        orbit = movingObject.getOrbit()
        src = orbit.getSrc()
        
        # Update the MovingObject table.
        db.setColumnLong('movingObjectId', movingObjectId)
        db.setColumnString('mopsStatus', movingObject.getStatus())
        h_v = movingObject.getH_v()
        if(h_v == None):
            db.setColumnToNull('h_v')
        else:
            db.setColumnDouble('h_v', h_v)
        g = movingObject.getG()
        if(g == None):
            db.setColumnToNull('g')
        else:
            db.setColumnDouble('g', g)
        db.setColumnDouble('q', orbit.getQ())
        db.setColumnDouble('e', orbit.getE())
        db.setColumnDouble('i', orbit.getI())
        db.setColumnDouble('node', orbit.getNode())
        db.setColumnDouble('argPeri', orbit.getArgPeri())
        db.setColumnDouble('timePeri', orbit.getTimePeri())
        for i in range(1, 22, 1):
            db.setColumnDouble('src%02d' %(i), src[i-1])
        db.insertRow()
    
    # Now add the MovingObject -> Tracklet info.
    # TODO: is it better to use 3 cursors and do all of this in the loop above?
    if(updateTrackletStatus):
        db.setTableForInsert('mops_MovingObjectToTracklet')
        for movingObject in movingObjects:
            movingObjectId = movingObject.getMovingObjectId()
            tracklets = movingObject.getTracklets()
            
            for tracklet in tracklets:
                db.setColumnLong('movingObjectId', movingObjectId)
                db.setColumnLong('trackletId', tracklet.getTrackletId())
                db.insertRow()
        
            # Finally, update the tracklet status. We assume that the caller has
            # set the status appropriately, of course.
            TrackletList.updateStatus(dbLocStr, tracklets)
    else:
        db.setTableForInsert('mops_MovingObjectToTracklet')
        for movingObject in movingObjects:
            movingObjectId = movingObject.getMovingObjectId()
            tracklets = movingObject.getTracklets()
            
            for tracklet in tracklets:
                db.setColumnLong('movingObjectId', movingObjectId)
                db.setColumnLong('trackletId', tracklet.getTrackletId())
                db.insertRow()
    return


def _getNextMovingObjectId(dbLocStr):
    # Since movingObjectId is autoincrement, it cannot be 0. It will get
    # incremented by 1 at the end of the function call...
    movingObjectId = 0
            
    # Connect to the database.
    db = SafeDbStorage()
    db.setRetrieveLocation(persistence.LogicalLocation(dbLocStr))
    
    db.setTableForQuery('MovingObject')
    db.outColumn('max(movingObjectId)', True)      # isExpr=True
    
    db.query()
    if(db.next() and not db.columnIsNull(0)):
        movingObjectId = db.getColumnByPosLong(0)
    db.finishQuery()
    return(movingObjectId + 1)





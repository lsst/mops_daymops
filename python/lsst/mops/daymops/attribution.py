"""
Generic functions for attribution and precovery.

Not finished: just an initial skeleton.
"""
import lib
from Tracklet import Tracklet

import lsst.daf.persistence as persistence
from SafeDbStorage import SafeDbStorage

import auton


# Constants



def getTonightsNightNumber(dbLocStr, utOffset):
    """
    Given a database location string, figure out the night number for tonight's
    DiaSources.
    
    The offset from UT (utOffset) is defined as
        local time = UT time - utOffset
    
    @param dbLocStr: database conection string.
    @param utOffset: offset from localtime to UT (in hours).
    """
    # Send the query
    # sql = select DiaSource.firstExposureTime from DiaSource, DiaSourceIDTonight
    #       where DiaSource.diaSourceId=DiaSourceIDTonight.DiaSourceId limit 1
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    db.setTableListForQuery(('DiaSource', 'DiaSourceIDTonight'))
    db.outColumn('DiaSource.firstExposureTime')
    db.outColumn('DiaSource.diaSourceId')
    db.setQueryWhere('DiaSource.diaSourceId=DiaSourceIDTonight.DiaSourceId')
    db.query()
    
    if(db.next()):
        mjd = db.getColumnByPosDouble(0)
    else:
        # We might not have any DiaSources for tonight. In that case raise an
        # exception.
        raise(Exception('Unable to determine night number: ' + \
                        'DiaSourceIDTonight is empty.'))
    db.finishQuery()
    del(db)
    return(lib.mjdToNightNumber(mjd, utOffset))











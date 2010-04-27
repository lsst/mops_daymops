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
    DIASources.
    
    The offset from UT (utOffset) is defined as
        local time = UT time - utOffset
    
    @param dbLocStr: database conection string.
    @param utOffset: offset from localtime to UT (in hours).
    """
    # Send the query
    # sql = select DIASource.taiMidPoint from DIASource, DIASourceIDTonight
    #       where DIASource.diaSourceId=DIASourceIDTonight.DIASourceId limit 1
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    db.setTableListForQuery(('DIASource', 'DIASourceIDTonight'))
    db.outColumn('DIASource.taiMidPoint')
    db.outColumn('DIASource.diaSourceId')
    db.setQueryWhere('DIASource.diaSourceId=DIASourceIDTonight.DIASourceId')
    db.query()
    
    if(db.next()):
        mjd = db.getColumnByPosDouble(0)
    else:
        # We might not have any DIASources for tonight. In that case raise an
        # exception.
        raise(Exception('Unable to determine night number: ' + \
                        'DIASourceIDTonight is empty.'))
    db.finishQuery()
    del(db)
    return(lib.mjdToNightNumber(mjd, utOffset))











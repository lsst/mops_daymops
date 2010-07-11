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











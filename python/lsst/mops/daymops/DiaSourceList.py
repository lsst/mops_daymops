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
Helper functions to create/retrieve a list of DIASource instances.
"""
from DayMOPSObject import DayMOPSObject
from DiaSource import DiaSource
import lib

import lsst.daf.persistence as persistence
from SafeDbStorage import SafeDbStorage




def diaSourceListForTonight(dbLocStr, sliceId=None, numSlices=None):
    """
    Use  sliceId and numSlices to implement some form of parallelism.
    
    @param dbLocStr: database connection string.
    @param sliceId: Id of the current Slice.
    @param numSlices: number of available slices (i.e. MPI universe size - 1)
    
    Return 
    Interator to the list of DIASource instances.
    """
    # Send the query.
    # sql: select d.diaSourceId, d.ra, d.decl, d.filterId, d.taiMidPoint, \
    #      d.obsCode, d.apFlux, d.apFluxErr, d.refMag from \
    #      DIASource d, DIASourceIDTonight t where \
    #      t.DIASourceId=d.diaSourceId;
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    db.setTableListForQuery(('DIASource', 'DIASourceIDTonight'))
    db.outColumn('DIASource.diaSourceId')
    db.outColumn('DIASource.ra')
    db.outColumn('DIASource.decl')
    db.outColumn('DIASource.filterId')
    db.outColumn('DIASource.taiMidPoint')
    db.outColumn('DIASource.obsCode')
    db.outColumn('DIASource.apFlux')
    db.outColumn('DIASource.apFluxErr')
    db.outColumn('DIASource.refMag')
    where = 'DIASource.diaSourceId=DIASourceIDTonight.DIASourceId'
    if(sliceId != None and numSlices > 1):
        where += ' and DIASource.diaSourceId %% %d = %d' \
                 %(numSlices, sliceId)
    db.setQueryWhere(where)
    db.query()
    
    # Fetch the results.
    # FIXME: Update these every time afw updates.
    while(db.next()):
        d = DiaSource()
        d.setDiaSourceId(db.getColumnByPosLong(0))
        d.setRa(db.getColumnByPosDouble(1))
        d.setDec(db.getColumnByPosDouble(2))
        d.setFilterId(db.getColumnByPosInt(3))
        d.setTaiMidPoint(float(db.getColumnByPosDouble(4)))
        d.setObsCode(db.getColumnByPosString(5))
        d.setApFlux(db.getColumnByPosDouble(6))
        d.setApFluxErr(db.getColumnByPosDouble(7))
        d.setRefMag(db.getColumnByPosDouble(8))
        yield(d)
    db.finishQuery()
    del(db)
    # return


# Methods used to compute some statistics from our DiaSource objects.
def computeVelocityStats(diaSources):
    """
    Compute basic velocity information based on the details of members of
    diaSources.
    
    @param diaSources: list of DIASource instances.
    
    Return 
    Velocity components and modulus, in deg/day computed from the input list of 
    DIASource positions and MJDs.
        [velRa, velDec, velModulus]
    """
    # Sort diaSources by MJD. Compute stats on first and last DiaSource. We 
    # assume linear motion from the first DiaSource to the last only for now.
    if(not diaSources or len(diaSources) < 2):
        return((None, None, None))
    
    diaSources.sort()
    first = diaSources[0]
    last = diaSources[-1]
    
    # Compute time distance in days.
    timeDistance = last.getTaiMidPoint() - first.getTaiMidPoint()
    if(not timeDistance):
        # FIXME: is this a good idea? Should we just return [0., 0., 0.]?
        raise(Exception('No temporal spread in DiaSources!'))
    
    # Compute spherical distance.
    distance = lib.sphericalDistance((first.getRa(), first.getDec()),
                                     (last.getRa(), last.getDec()))
    return([d / timeDistance for d in distance])

def getTimeSpan(diaSources):
    """
    Return the time span in day of diaSources. The time span is 
    defined as
        max(d.getTaiMidPoint())-min(d.getTaiMidPoint) for d in diaSources
    
    @param diaSources: list of DIASource instances.
    
    Return
    The time spanned (in days) by the input DIASource instances.
    """
    if(not diaSources):
        return(None)
    
    times = [d.getTaiMidPoint() for d in diaSources]
    times.sort()
    return(times[-1] - times[0])






























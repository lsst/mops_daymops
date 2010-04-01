"""
Helper functions to create/retrieve a list of DiaSource instances.
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
    Interator to the list of DiaSource instances.
    """
    # Send the query.
    # sql: select d.diaSourceId, d.ra, d.decl, d.filterId, d.exposureStartTime, \
    #      d.obsCode, d.apFlux, d.apFluxSigma, d.psfFlux from \
    #      DiaSource d, DiaSourceIDTonight t where \
    #      t.DiaSourceId=d.diaSourceId;
    db = SafeDbStorage()
    db.setPersistLocation(persistence.LogicalLocation(dbLocStr))
    db.setTableListForQuery(('DiaSource', 'DiaSourceIDTonight'))
    db.outColumn('DiaSource.diaSourceId')
    db.outColumn('DiaSource.ra')
    db.outColumn('DiaSource.decl')
    db.outColumn('DiaSource.filterId')
    #jmyers: it looks like DiaSource has taiMidPoint even in newest AFW
    # so just cheat and store this there... TBD: find out what's
    # going on with AFW/schema inconsistency and straighten them out!
    db.outColumn('DiaSource.exposureStartTime')
    db.outColumn('DiaSource.apFlux')
    db.outColumn('DiaSource.apFluxSigma')
    #jmyers: another weird inconsistency between schema and AFW; AFW
    # uses refmag and schema does not.  It doesn't really matter
    # since we don't use magnitude data in MOPS.  TBD: Same as above.
    db.outColumn('DiaSource.psfFlux')
    where = 'DiaSource.diaSourceId=DiaSourceIDTonight.DiaSourceId'
    if(sliceId != None and numSlices > 1):
        where += ' and DiaSource.diaSourceId %% %d = %d' \
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
        # TBD: change this!
        d.setTaiMidPoint(float(db.getColumnByPosDouble(4)))
        d.setApFlux(db.getColumnByPosDouble(6))
        d.setApFluxErr(db.getColumnByPosDouble(7))
        # TBD: change this too!
        d.setPsfFlux(db.getColumnByPosDouble(8))
        yield(d)
    db.finishQuery()
    del(db)
    # return


# Methods used to compute some statistics from our DiaSource objects.
def computeVelocityStats(diaSources):
    """
    Compute basic velocity information based on the details of members of
    diaSources.
    
    @param diaSources: list of DiaSource instances.
    
    Return 
    Velocity components and modulus, in deg/day computed from the input list of 
    DiaSource positions and MJDs.
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
    timeDistance = last.getExposureStartTime() - first.getExposureStartTime()
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
    
    @param diaSources: list of DiaSource instances.
    
    Return
    The time spanned (in days) by the input DiaSource instances.
    """
    if(not diaSources):
        return(None)
    
    times = [d.getTaiMidPoint() for d in diaSources]
    times.sort()
    return(times[-1] - times[0])






























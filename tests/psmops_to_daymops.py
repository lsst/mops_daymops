#!/usr/bin/env python
"""
Take the content of a PS-MOPS simulation database and massage it into a 
LSST-MOPS initial database. This assumes that the both databases exist and that
the LSST-MOPS database is already setup but empty.


Usage
    shell> psmops_to_daymops.py <psmops_db> <lsst_db>

Example
    shell> ./psmops_to_daymops.py psmops_onelunation daymops_onelunation


Revision History
April 23, 2009: First version (F. Pierfederici)
"""
import math
import sys
import MySQLdb as DBI
from lsst.daf.base import DateTime



# Constants
DB_HOST = 'localhost'
DB_USER = 'jmyers'
DB_PASSWORD = 'jmyers'

FOV_RADIUS = 3.5 / 2.                       # Size of the FoV radius in degrees.

DEG_TO_RAD = math.pi / 180.
RAD_TO_DEG = 1. / DEG_TO_RAD

USAGE = '''Usage
    shell> psmops_to_daymops.py <psmops_db> <lsst_db>

Example
    shell> ./psmops_to_daymops.py psmops_onelunation daymops_onelunation

'''



def _toSqlString(x):
    """
    Convert to string taking care of the None->NULL conversion.
    """
    if(x == None):
        return('NULL')
    return(str(x))

    
def _populateTable(readCur, readSql, writeCur, writeSql):
    print "fetching records..."
    numRes = readCur.execute(readSql)
    print "done."
    row = readCur.fetchone()
    i = 0
    while(row):
        args = tuple([_toSqlString(x) for x in row])
        n = writeCur.execute(writeSql % args)
        i += 1
        row = readCur.fetchone()
    return(i)


def populateSSM(psmopsConn, daymopsConn):    
    psmopsCur = psmopsConn.cursor()
    daymopsCur = daymopsConn.cursor()
    
    # ssm_desc -> mops_SSMDesc
    # First clear the table.
    daymopsCur.execute('delete from `mops_SSMDesc`')
    daymopsConn.commit()
    
    # Convert the data.
    readSql = 'select `desc_id`, `prefix`, `description` from `ssm_desc`'
    writeSql = 'insert into `mops_SSMDesc` \
                (`ssmDescId`, `prefix`, `description`) values (%s, "%s", "%s")'
    print "writing to mops_SSMDesc"
    n = _populateTable(psmopsCur, readSql, daymopsCur, writeSql)
    daymopsConn.commit()
    
    
    # ssm -> mops_SSM
    # First clear the table.
    daymopsCur.execute('delete from `mops_SSM`')
    daymopsConn.commit()
    
    # Convert the data.
    readSql = 'select `ssm_id`, `q`, `e`, `i`, `node`, \
               `arg_peri`, `time_peri`, `epoch`, `h_v`, `h_ss`, `g`, `albedo`, \
               `object_name`, `desc_id` from `ssm`'
    writeSql = 'insert into `mops_SSM` (`ssmId`, `q`, `e`, `i`, `node`, \
                `argPeri`, `timePeri`, `epoch`, `h_v`, `h_ss`, `g`, `albedo`, \
                `ssmObjectName`,`ssmDescId`) values (%s, %s, %s, %s, %s, \
                %s, %s, %s, %s, %s, %s, %s, \
                "%s", %s)'
    print "writing to mops_SSM"
    n = _populateTable(psmopsCur, readSql, daymopsCur, writeSql)
    daymopsConn.commit()
    return


def populateExposure(psmopsConn, daymopsConn):
    psmopsCur = psmopsConn.cursor()
    daymopsCur = daymopsConn.cursor()
    
    # filter_info -> prv_Filter
    # Forget about it for now since it would require populate the focal plane
    # table as well... we will just assume all PS-MOPS fields have are in one 
    # filter and assign it aa arbitrary ID of 0.
    # First clear the table.
    filterId = 0
    
    # `fields` -> Raw_Amp_Exposure + Raw_CCD_Exposure + Raw_FPA_Exposure
    # First clear the table.
    daymopsCur.execute('delete from `Raw_Amp_Exposure`')
    daymopsConn.commit()
    daymopsCur.execute('delete from `Raw_CCD_Exposure`')
    daymopsConn.commit()
    daymopsCur.execute('delete from `Raw_FPA_Exposure`')
    daymopsConn.commit()
    
    # Convert the data.
    readSql = 'select `field_id`, `epoch_mjd`, `ra_deg`, `dec_deg`, \
               `time_start`, `time_stop`, `limiting_mag`, \
               `ra_sigma`, `dec_sigma`, `obscode` from `fields`'
    writeFPASql = 'insert into Raw_FPA_Exposure (`rawFPAExposureId`, \
                   `filterId`, `ra`, `decl`, `obsDate`, `tai`, `ra_ll`, \
                   `dec_ll`, `ra_lr`, `dec_lr`, `ra_ul`, `dec_ul`, \
                   `ra_ur`, `dec_ur`) values (%s, \
                   %s, %s, %s, "%s", %s, %s, \
                   %s, %s, %s, %s, %s, \
                   %s, %s)'
    writeCCDSql = 'insert into Raw_Amp_Exposure (rawAmpExposureId, \
                   rawCCDExposureId, taiObs, expTime) values (%s, \
                   %s, "%s", %s)'
    writeAmpSql = 'insert into Raw_CCD_Exposure (rawCCDExposureId, \
                   rawFPAExposureId, filterId, ra, decl, dateObs, taiObs, \
                   mjdObs, expTime) values (%s, \
                   %s, %s, %s, %s, "%s", "%s", \
                   %s, %s)'
                
    # Massage the data.
    numRes = psmopsCur.execute(readSql)
    row = psmopsCur.fetchone()
    while(row):
        (field_id, epoch_mjd, ra, dec, time_start, time_stop, 
         limiting_mag, ra_sigma, dec_sigma, obscode) = row
        
        expTime = abs(time_stop - time_start)
        
        # Convert the UTC MJDs to TAI nsecs and TAI datetimes.
        t = DateTime(epoch_mjd, DateTime.UTC)
        tai = t.nsecs(DateTime.TAI)
        mjdTai = t.mjd(DateTime.TAI)
        # Some tables have taiObs, some have dateObs. Since we do not use those,
        # we assume that they are the same... which is not true (they differ by
        # a few seconds).
        obsDate = t.toString()
        # Fix it so that MySQL is happy!
        obsDate = obsDate.replace('T', ' ')
        obsDate = obsDate.replace('Z', '')
        obsDate = obsDate[:19]
        
        decRadians = dec * DEG_TO_RAD
        dec_ll = dec - FOV_RADIUS
        dec_lr = dec - FOV_RADIUS
        dec_ul = dec + FOV_RADIUS
        dec_ur = dec + FOV_RADIUS
        ra_ll = ra - FOV_RADIUS * math.cos(dec_ll)
        ra_lr = ra + FOV_RADIUS * math.cos(dec_lr)
        ra_ul = ra - FOV_RADIUS * math.cos(dec_ul)
        ra_ur = ra + FOV_RADIUS * math.cos(dec_ur)
        
        # Write FPA.
        n = daymopsCur.execute(writeFPASql % (str(field_id),
                                              str(filterId),
                                              str(ra),
                                              str(dec),
                                              obsDate,
                                              str(tai),
                                              str(ra_ll),
                                              str(dec_ll),
                                              str(ra_lr),
                                              str(dec_lr),
                                              str(ra_ul),
                                              str(dec_ul),
                                              str(ra_ur),
                                              str(dec_ur)))
        # Write CCD.
        n = daymopsCur.execute(writeCCDSql % (str(field_id),
                                              str(field_id),
                                              obsDate,
                                              str(expTime)))
        # Write Amp.
        n = daymopsCur.execute(writeAmpSql % (str(field_id),
                                              str(field_id),
                                              str(filterId),
                                              str(ra),
                                              str(dec),
                                              obsDate,
                                              obsDate,
                                              str(mjdTai),
                                              str(expTime)))
        row = psmopsCur.fetchone()
    daymopsConn.commit()
    return


def populateSources(psmopsConn, daymopsConn):
    psmopsCur = psmopsConn.cursor()
    daymopsCur = daymopsConn.cursor()

    # `detections` -> DIASource
    # First clear the table.
    daymopsCur.execute('delete from `DIASource`')
    daymopsConn.commit()
    
    # Convert the data.
    readSql = 'select `det_id`, `field_id`, `ra_deg`, `dec_deg`, `epoch_mjd`, \
               `mag`, `ref_mag`, `filter_id`, `s2n`, `ra_sigma_deg`, \
               `dec_sigma_deg`, `mag_sigma`, `orient_deg`, `length_deg`, \
               `ssm_id`, `obscode`, `is_synthetic` from detections, ssm where \
               ssm.object_name=detections.object_name'
    writeDIASql = 'insert into DIASource (diaSourceId, ampExposureId, \
                   filterId, ssmId, ra, raErrForDetection, decl, \
                   declErrForDetection, taiMidPoint, apFlux, apFluxErr, \
                   refMag, obsCode, isSynthetic) values (%s, %s, \
                   %s, %s, %s, %s, %s, \
                   %s, %s, %s, %s, \
                   %s, "%s", %s)'
    
    # Massage the data.
    numRes = psmopsCur.execute(readSql)
    row = psmopsCur.fetchone()
    while(row):
        (detId, fieldId, ra, dec, epochMjd, mag, refMag, filterId, s2n,
         raErr, decErr, magErr, orient, length, ssmId, obsCode, 
         isSynthetic) = row
        
        # Again filter_id is 0 for now.
        filterId = 0
        
        # Convert the UTC MJDs to TAI nsecs and TAI datetimes.
        print "source has epochMJD = ", epochMJD
        t = DateTime(epochMjd, DateTime.UTC)
        taiMidPoint = t.nsecs(DateTime.TAI)
        print "calculated associated taiMidPoint = ", taiMidPoint
        
        # Compute flux and flux error.
        # m = m0 -2.5Log(f)
        # f = 10^-0.4(m-m0)
        # fErr = magErr * ln(10) * 10^-0.4(m-m0)
        flux = 10**(-0.4 * (mag - refMag))
        fluxErr = flux * magErr * math.log(10)
        
        # Write FPA.
        n = daymopsCur.execute(writeDIASql % (str(detId),
                                              str(fieldId),
                                              str(filterId),
                                              str(ssmId),
                                              str(ra),
                                              str(raErr),
                                              str(dec),
                                              str(decErr),
                                              str(taiMidPoint),
                                              str(flux),
                                              str(fluxErr),
                                              str(refMag),
                                              obsCode,
                                              isSynthetic))
        row = psmopsCur.fetchone()
    daymopsConn.commit()
    return




if(__name__ == '__main__'):
    # Get the database names.
    if(len(sys.argv) != 3):
        sys.stderr.write(USAGE)
        sys.exit(1)
    psmopsDb = sys.argv[1]
    daymopsDb = sys.argv[2]
    
    
    # Open two connections to the DB.
    psmopsConn = DBI.connect(host=DB_HOST, user=DB_USER, passwd=DB_PASSWORD,
                             db=psmopsDb)
    daymopsConn = DBI.connect(host=DB_HOST, user=DB_USER, passwd=DB_PASSWORD,
                             db=daymopsDb)
    
    
    # SSM
    populateSSM(psmopsConn, daymopsConn)
    
    # Exposure tables.
    populateExposure(psmopsConn, daymopsConn)
    
    # detections -> DIASource
    populateSources(psmopsConn, daymopsConn)
    
    
    
    

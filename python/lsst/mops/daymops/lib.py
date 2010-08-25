"""
Utility functions for DayMOPS

Here you will find generic funtions used throughout DayMOPS.
"""
from lsst.daf.base import DateTime

import math



# Constants
DEFAULT_UT_OFFSET = -4. / 24.
DEG_TO_RAD = math.pi / 180.
RAD_TO_DEG = 1. / DEG_TO_RAD
RAD_TO_DEG_SQ = RAD_TO_DEG**2.
LOG10 = math.log(10)
K = -2.5 / LOG10
pi2 = math.pi / 2.




def nsecToNightNumber(nsecs, utOffset=DEFAULT_UT_OFFSET):
    """
    Given a date in UT TAI nsecs and an offset (in hours) from UT, compute the
    corresponding night number.
    
    Night numbers are defined as
        int(DateTime(nsecs).mjd - .5 + utOffset / 24.)
    which means that the night number increases by 1 at the local noon. All 
    nsecs between two consecutive local noons correspond to the same night 
    number.
    
    The offset from UT (utOffset) is defined as
        local time = UT time - utOffset
    
    @param nsecs: the date to convert (in UT TAI nsecs) to night number.
    @param utOffset: offset from localtime to UT (in hours).
    """
    return(mjdToNightNumber(DateTime(nsecs).mjd(DateTime.TAI), utOffset))

    
def nightNumberToNsecRange(nightNumber, utOffset=DEFAULT_UT_OFFSET):
    """
    Given a night number and an offset (in hours) from UT, compute the 
    corresponding UT TAI nsecs range in the form (nsecMin, nsecMax).
    
    Night numbers are defined as
        int(DateTime(nsecs).mjd - 0.5 + utOffset / 24.0)
    which means that the night number increases by 1 at the local noon. All 
    nsecs between two consecutive local noons correspond to the same night 
    number.
    
    The offset from UT (utOffset) is defined as
        local time = UT time - utOffset
    
    @param nightNumber: the night number to convert to UT TAI nsecs.
    @param utOffset: offset from localtime to UT (in hours).
    """
    return([DateTime(t).nsecs(DateTime.TAI) \
            for t in nightNumberToMjdRange(nightNumber, utOffset)])


def mjdToNightNumber(mjd, utOffset=DEFAULT_UT_OFFSET):
    """
    Given a date in UT TAI MJD and an offset (in hours) from UT, compute the
    corresponding night number.
    
    Night numbers are defined as
        int(mjd - .5 + utOffset / 24.)
    which means that the night number increases by 1 at the local noon. All MJDs
    between two consecutive local noons correspond to the same night number.
    
    The offset from UT (utOffset) is defined as
        local time = UT time - utOffset
    
    @param mjd: the date to convert (in UT TAI MJD) to night number.
    @param utOffset: offset from localtime to UT (in hours).
    """
    return(int(mjd - .5 + utOffset / 24.))

    
def nightNumberToMjdRange(nightNumber, utOffset=DEFAULT_UT_OFFSET):
    """
    Given a night number and an offset (in hours) from UT, compute the 
    corresponding UT TAI MJD range in the form (mjdMin, mjdMax).
    
    Night numbers are defined as
        int(mjd - 0.5 + utOffset / 24.0)
    which means that the night number increases by 1 at the local noon. All MJDs
    between two consecutive local noons correspond to the same night number.
    
    The offset from UT (utOffset) is defined as
        local time = UT time - utOffset
    
    @param nightNumber: the night number to convert to UT TAI MJD.
    @param utOffset: offset from localtime to UT (in hours).
    """
    offset = utOffset / 24.
    mjdMin = nightNumber + .50000 - offset
    mjdMax = nightNumber + 1.49999 - offset
    return(mjdMin, mjdMax)


def fluxToMag(flux, fluxErr, refMag=30.):
    """
    Given a flux, its error and a reference magnitude, transform fluxes in mags 
    and flux errors in mag errors.
    
    We use these formulas:
        mag = refMag - 2.5 * Log(flux)
        magErr = (-2.5 / ln(10)) * (fluxErr / flux)
    
    @param flux: flux to convert to mag
    @param fluxErr: error on flux
    @param refMag: reference magnitude (usually 30.)
    
    Return
        (mag, magErr)
    """
    mag = refMag - 2.5 * math.log(flux, 10)
    magErr = K * fluxErr / flux
    return(mag, magErr)


def magToFlux(mag, magErr, refMag=30.):
    """
    Given a mag, its error and a reference magnitude, transform mags in fluxes 
    and mag errors in flux errors.
    
    We use these formulas:
        flux = 10**(-0.4 * (mag - refMag))
        fluxErr = ln(10) * magErr * flux
    
    @param mag: mag to convert to flux
    @param magErr: error on mag
    @param refMag: reference magnitude (usually 30.)
    
    Return
        (flux, fluxErr)
    """
    flux = 10.**(-0.4 * (mag - refMag))
    fluxErr = flux * magErr * LOG10
    return(flux, fluxErr)


def sphericalDistance(point1, point2):
    """
    Compute the angular distance from point1 to point2 on the sky. All angles 
    are in degrees.
    
    @param point1: (RA, Dec) of the first point (degrees).
    @param point2: (RA, Dec) of the second point (degrees).
    
    Return (distanceRa, distanceDec, distanceModulus) (degrees)
    """
    if(abs(point1[1] - point2[1]) < 1.):
        return(_fastSphericalDistance(point1, point2))
    return(_slowSphericalDistance(point1, point2))


def _fastSphericalDistance(point1, point2):
    """
    Only valid for small distances.
    """
    # Convert everything to radians.
    p1 = [x * DEG_TO_RAD for x in point1]
    p2 = [x * DEG_TO_RAD for x in point2]
    
    # Compute the distances.
    cosDec = math.cos((p1[1] + p2[1]) / 2.)
    
    raDist = (p2[0] - p1[0]) * cosDec * RAD_TO_DEG
    decDist = (p2[1] - p1[1]) * RAD_TO_DEG
    totDist = math.sqrt(raDist**2 + decDist**2)
    return((raDist, decDist, totDist))


def _slowSphericalDistance(point1, point2):
    """
    General formula.
    """
    [ra1, dec1] = [x * DEG_TO_RAD for x in point1]
    [ra2, dec2] = [x * DEG_TO_RAD for x in point2]
    h = (pi2 - dec1)
    k = (pi2 - dec2)
    w = math.cos(ra2-ra1)
    cosh = math.cos(h)
    cosk = math.cos(k)
    sinh = math.sin(h)
    sink = math.sin(k)
    cosh_cosk = cosh * cosk
    sinh_sink = sinh * sink
    
    raDist = math.acos(cosh**2 + w * sinh**2) * RAD_TO_DEG_SQ * math.cos(dec1)
    while(raDist > 180.):
        raDist -= 180.
    decDist = math.acos(cosh_cosk + sinh_sink) * RAD_TO_DEG
    while(decDist > 180.):
        decDist -= 180.
    totDist = math.acos(cosh_cosk + sinh_sink * w) * RAD_TO_DEG
    return((raDist, decDist, totDist))














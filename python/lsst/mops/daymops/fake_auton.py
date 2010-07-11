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
This is just for testing: do not use in production!
"""
import random
import warnings
warnings.warn('Please use the real Auton module')





def findtracklets(detections, *args, **kwargs):
    """
    Just return random associations of detections as a list of the form
        [[detId1, detId2, ...], ]
    
    Input detections are a list of 
        [id, mjd, ra, dec, mag, obscode, objName, trailLength, trailAngle]
    
    @param detections: list of input detections.
    All other parameters are ignored.
    """
    if(len(detections) < 2):
        return([])
    
    diaSourceIds = [d[0] for d in detections]
    diaSourceMjds = [d[1] for d in detections]
    
    # Choose a random number of groups to return, between 0 and 
    # len(diaSourceIds) / 2 (2 detections per tracklet).
    n = int(random.random() * len(diaSourceIds) / 2.)
    
    # Now create n groups of indeces, each 2 to 5 elements long.
    # Also make sure that no detections with the same mjds end up in the same 
    # group.
    groups = []
    for i in range(n):
        tot = int(random.uniform(2, 5))
        indeces = []
        mjds = []
        for x in range(tot):
            idx = int(random.random()*(len(detections)-1))
            mjd = diaSourceMjds[idx]
            while(idx in indeces or mjd in mjds):
                idx = int(random.random()*(len(detections)-1))
                mjd = diaSourceMjds[idx]
            indeces.append(idx)
            mjds.append(mjd)
        groups.append([diaSourceIds[idx] for idx in indeces])
    return(groups)


def linktracklets(detections, *args, **kwargs):
    """
    Just return random associations of tracklets as a list of the form
        [[trackletId1, trackletId2, ...], ]
    
    Input detections list havs the form
        [(trackletId, detMjd, detRa, detDec, detMag, detObsCode, trackletName),]
    Several entries share the same tracletId, since each tracklet has >= 2 
    DiaSources.
    
    @param detections: list of input tracklets/detections.
    All other parameters are ignored.
    """
    # Get the number of unique trackletIds, this is the only thing we care 
    # about.
    trackletIds = []
    for det in detections:
        _id = det[0]
        if(_id not in trackletIds):
            trackletIds.append(_id)
    
    # We want at least 3 tracklets per track.
    if(len(trackletIds) < 3):
        return([])
    
    # Choose a random number of groups to return, between 0 and 
    # len(trackletIds) / 3 (3 detections per tracklet).
    n = int(random.random() * len(trackletIds) / 3.)
    
    # Now create n groups of indeces, each 3 to 10 elements long.
    groups = []
    for i in range(n):
        tot = int(random.uniform(3, 10))
        indeces = []
        for x in range(tot):
            idx = int(random.random()*(len(trackletIds)-1))
            while(idx in indeces):
                idx = int(random.random()*(len(trackletIds)-1))
            indeces.append(idx)
        groups.append([trackletIds[idx] for idx in indeces])
    return(groups)











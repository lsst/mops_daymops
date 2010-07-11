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
Class to represent a Tracklet object.
"""
from DayMOPSObject import DayMOPSObject
import DiaSourceList
import lsst.daf.persistence as persistence


# Constants
STATUS = {'UNATTRIBUTED':   'U',
          'ATTRIBUTED':     'A',
          'KILLED':         'K'}


class Tracklet(DayMOPSObject):
    def __init__(self, 
                 trackletId=None, 
                 status=STATUS['UNATTRIBUTED'],
                 diaSources=[]):
        self._trackletId = trackletId
        self._status = status
        self.setDiaSources(diaSources)
        return
    
    def __str__(self):
        return('Tracklet(trackletId=%d)' %(self._trackletId))
    
    def setDiaSources(self, diaSources):
        """
        Set the internal diaSources.
        
        @param diaSources: list of DiaSource instances.
        """
        self._diaSources = diaSources
        
        # Update velocity data.
        if(self._diaSources):
            (self._velRa, 
             self._velDec, 
             self._velTot) = DiaSourceList.computeVelocityStats(self._diaSources)
        return
    
    # Aliases
    def setVelDecl(self, v):
        return(self.setVelDec(v))
    
    
        
        
























        
        
    






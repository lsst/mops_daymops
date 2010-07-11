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
Class to represent a MovingObject object.
"""
from DayMOPSObject import DayMOPSObject
import TrackletList



# Constants
STATUS = {'NEW':            'N',
          'MERGED':         'M',
          'IOD FAILED':     'I',
          'DIFF FAILED':    'F',
          'OK':             'Y',
          'PRELIMINARY':    'T'}


class MovingObject(DayMOPSObject):
    """
    Representation of a MovingObject.
    
    This, together with the Orbit class, fully describe the subset of the
    MovingObject table of interset to MOPS.
    """
    def __init__(self, 
                 movingObjectId=None, 
                 status=STATUS['NEW'],
                 orbit=None,
                 h_v=None,
                 g=.150,
                 tracklets=[],
                 arcLength=None):
        super(MovingObject, self).__init__()
        
        self._movingObjectId = movingObjectId
        self._status = status
        self._orbit = orbit
        self._h_v = h_v
        self._g = None
        self.setG(g)
        self._arcLength = arcLength
        
        # This updates the arcLength...
        self._tracklets = None
        self.setTracklets(tracklets)
        return
    
    def setG(self, g=.150):
        """
        Set self._g to the desired value. If g=None, then use the default value
        of .150
        """
        if(g == None):
            self._g = .150
        else:
            self._g = g
        return
    
    def setTracklets(self, tracklets):
        """
        Autmatically compute the arc length, updating it.
        """
        self._tracklets = tracklets
        if(tracklets):
            self._arcLength = TrackletList.getArcLength(tracklets)
        return
    
    def __str__(self):
        return('MovingObject(movingObjectId=%d)' %(self._movingObjectId))
    
    # Aliases
    def setArcLengthDays(self, length):
        """
        Alias to self.setArcLength()
        """
        return(self.setArcLength(length))
    
    def getArcLengthDays(self):
        """
        Alias to self.getArcLength()
        """
        return(self.getArcLength())
    
        
    
    
    
    
        
        
























        
        
    






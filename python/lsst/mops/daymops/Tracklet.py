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
    
    
        
        
























        
        
    






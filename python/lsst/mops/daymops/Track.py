"""
Class to represent a Track object.

Note that Track ojects do not really exist in the database. What one finds in 
the database are just groupings of trackletIds. Each group corresponds to a 
Track "instance".
"""
from DayMOPSObject import DayMOPSObject
import lsst.daf.persistence as persistence


# Constants
# None yet.


class Track(DayMOPSObject):
    def __init__(self, 
                 trackId=None, 
                 tracklets=[]):
        self._trackId = trackId
        self._tracklets = tracklets
        return
    
    def __str__(self):
        return('Track(trackId=%d)' %(self._trackId))

    
    
        
        
























        
        
    






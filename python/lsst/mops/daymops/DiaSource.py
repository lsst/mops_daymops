"""
Class to represent a DIASource object.

This is monkeypatching the corresponding C++ class.
"""
import lsst.afw.detection as detection




# This should realy left-inherit from DayMOPSObject, but that doesn't work with
# SWIG :-(
class DiaSource(detection.DiaSource):
    """
    Opaque layer above the DIASource and DIASOurceIDTonight tables.
    
    We compare DiaSources by their MJD alone, not their spatial location.
    """
    _refMag = None
    # FIXME: to be obsoleted by fix to bug #796
    _obsCode = ''
    
    def getRefMag(self):
        """
        Getter for refMag
        """
        return(self._refMag)
    
    def setRefMag(self, mag):
        """
        Setter for refMag
        """
        self._refMag = mag
        return
    
    def getObsCode(self):
        """
        Getter for obsCode
        """
        return(self._obsCode)
    
    def setObsCode(self, obsCode):
        """
        Setter for obsCode
        """
        if(not isinstance(obsCode, str) or len(obsCode) != 3):
            raise(SyntaxError('obsCode has to be a 3-letter string.'))
        self._obsCode = str(obsCode)
        return
    
    # Aliases
    def setDecl(self, dec):
        """
        Alias for self.setDec()
        """
        return(self.setDec(dec))
    
    def getDecl(self):
        """
        Alias for self.getDec()
        """
        return(self.getDec())
        
    # Comparison by MJD.
    def __lt__(self, other):
        if(other == None):
            return(False)
        return(self.getTaiMidPoint() < other.getTaiMidPoint())
    
    def __le__(self, other):
        if(other == None):
            return(False)
        return(self.getTaiMidPoint() <= other.getTaiMidPoint())
    
    def __eq__(self, other):
        if(other == None):
            return(False)
        return(self.getTaiMidPoint() == other.getTaiMidPoint())
    
    def __ne__(self, other):
        if(other == None):
            return(True)
        return(self.getTaiMidPoint() != other.getTaiMidPoint())

    def __gt__(self, other):
        if(other == None):
            return(False)
        return(self.getTaiMidPoint() > other.getTaiMidPoint())
    
    def __ge__(self, other):
        if(other == None):
            return(False)
        return(self.getTaiMidPoint() >= other.getTaiMidPoint())
    
    
    
    
    
    
    
    
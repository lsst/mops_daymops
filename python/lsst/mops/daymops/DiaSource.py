"""
Class to represent a DiaSource object.

This is monkeypatching the corresponding C++ class.
"""
import lsst.afw.detection as detection




# This should realy left-inherit from DayMOPSObject, but that doesn't work with
# SWIG :-(
class DiaSource(detection.DiaSource):
    """
    Opaque layer above the DiaSource and DIASOurceIDTonight tables.
    
    We compare DiaSources by their MJD alone, not their spatial location.
    """
    _psfFlux = None
    # FIXME: to be obsoleted by fix to bug #796
    _obsCode = ''
    
    def getPsfFlux(self):
        """
        Getter for psfFlux
        """
        return(self._psfFlux)
    
    def setPsfFlux(self, mag):
        """
        Setter for psfFlux
        """
        self._psfFlux = mag
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

    #alias by jmyers - TBD: get rid of this once AFW and Schema agree on which we're using!
    def getExposureStartTime(self):
        """ 
        alias for getTaiMidPoint() - one of these will become deprecated soon.
        """
        return self.getTaiMidPoint()
    def setExposureStartTime(self):
        """
        alias for setTaiMidPoint() - one of thse will become deprecated soon.
        """
        return self.setTaiMidPoint()
        
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
    
    
    
    
    
    
    
    

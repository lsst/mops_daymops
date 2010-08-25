"""
Since at the time of writing (July 23, 2009) the return value of the various
getColumnByPos*() methods in daf_persistence.DbStorage is undefined when the
coresponding value in the DB is NULL, we are monkeypatching the class to be more
safe.
"""
from lsst.daf.persistence import DbStorage



class SafeDbStorage(DbStorage):
    def _safeFetcher(self, dataType, columnIndex):
        if(self.columnIsNull(columnIndex)):
            return(None)
        m = getattr(super(SafeDbStorage, self), 'getColumnByPos%s' %(dataType))
        return(m(columnIndex))
    
    def getColumnByPosChar(self, i):
        return(self._safeFetcher('Char', i))
    
    def getColumnByPosShort(self, i):
        return(self._safeFetcher('Short', i))
    
    def getColumnByPosInt(self, i):
        return(self._safeFetcher('Int', i))
    
    def getColumnByPosLong(self, i):
        return(self._safeFetcher('Long', i))
        
    def getColumnByPosInt64_t(self, i):
        return(self._safeFetcher('Int64_t', i))
    
    def getColumnByPosFloat(self, i):
        return(self._safeFetcher('Float', i))
    
    def getColumnByPosDouble(self, i):
        return(self._safeFetcher('Double', i))
    
    def getColumnByPosString(self, i):
        return(self._safeFetcher('String', i))
    
    def getColumnByPosBool(self, i):
        return(self._safeFetcher('Bool', i))




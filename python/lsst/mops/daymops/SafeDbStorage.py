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




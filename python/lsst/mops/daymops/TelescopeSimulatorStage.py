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
Simulate the flow of data into DIASourceIDTonight

Used for testing purposes only.
"""
from DayMOPSStage import DayMOPSStage
import lib

import lsst.daf.persistence as persistence
from SafeDbStorage import SafeDbStorage





class TelescopeSimulatorStage(DayMOPSStage):
    """
    Simulate the telescope observing for one full night and populate the 
    DIASourceIDTonight database table with the appropriate DIASource IDs.
    
    This Stage should not be used in operations.
    """
    def __init__(self, stageId=-1, policy=None, startFromNightNumber=0):
        """
        Standard Stage initializer.
        """
        super(TelescopeSimulatorStage, self).__init__(stageId, policy)
        
        # Get database details from policy
        self.dbLocStr = self.getValueFromPolicy('database')
        
        # Times in .paf files are always in days, we need hours.
        self.utOffset = self.getValueFromPolicy('utOffset', 
                                                lib.DEFAULT_UT_OFFSET) * 24.
        
        
        
        # Fetch the oldest DIASource time (in UT TAI nsec).
        # Connect to the database.
        db = SafeDbStorage()
        db.setRetrieveLocation(persistence.LogicalLocation(self.dbLocStr))
        
        # If we do not have data between those two times, go to the next 
        # available times.
        # sql: select taiMidPoint from DIASource order by taiMidPoint limit 1;
        db.setTableForQuery('DIASource')
        db.outColumn('min(taiMidPoint)', True)      # isExpr=True
        db.outColumn('max(taiMidPoint)', True)      # isExpr=True
        db.query()
        if(db.next()):
            tMin = db.getColumnByPosDouble(0)
            tMax = db.getColumnByPosDouble(1)
            db.finishQuery()
        else:
            db.finishQuery()
            raise(Exception('No data in %s/DIASource' %(self.dbLocStr, )))
        del(db)
        
        # Now determine the starting night number.
        self.lastProcessedNight = max(int(startFromNightNumber),
                                      lib.mjdToNightNumber(tMin))
        self.lastProcessedNight -= 1
        self.theLatestPossibleNight = lib.mjdToNightNumber(tMax)
        self.logIt('INFO', 'Processing nights %d to %d' \
                   %(self.lastProcessedNight + 1, self.theLatestPossibleNight))
        return
    
    def preprocess(self):
        """
        Copy all the DIASources from night to night into the DIASourceIDTonight
        table. Start from the oldest night not grater than 
        self.lastProcessedNight.
        """
        # Call the superclass preprocess() method.
        super(TelescopeSimulatorStage, self).preprocess()
        
        # Update the night counter.
        self.lastProcessedNight += 1
        
        # Wipe the table clean.
        db = SafeDbStorage()
        db.setPersistLocation(persistence.LogicalLocation(self.dbLocStr))
        db.startTransaction()
        db.executeSql('delete from DIASourceIDTonight')
        db.endTransaction()
        del(db)
        
        # Do we have any more data to process?
        if(self.lastProcessedNight > self.theLatestPossibleNight):
            self.logIt('INFO', 'Closing up shop: nothing else to do!')
            self._shutdownPipeline()
            return
        
        self.logIt('INFO', 'Processing night %d.' %(self.lastProcessedNight))
        
        # Determine the nsec range from the night number.
        (tMin, tMax) = lib.nightNumberToMjdRange(self.lastProcessedNight,
                                                 self.utOffset)
        
        # Copy the corresponding DIASource IDs into DIASourceIDTonight.
        # sql: insert into DIASourceIDTonight (DIASourceId) \
        #      (select diaSourceId from DIASource \
        #       where taiMidPoint between %f and %f);
        sql = 'insert into DIASourceIDTonight (DIASourceId) \
               (select diaSourceId from DIASource \
                where taiMidPoint between %f and %f)'
        db = SafeDbStorage()
        db.setRetrieveLocation(persistence.LogicalLocation(self.dbLocStr))
        db.startTransaction()
        db.executeSql(sql %(tMin, tMax))
        db.endTransaction()
        del(db)
        self.logIt('INFO', 'Inserted DiaSource IDs into DIASourceIDTonight')
        self.outputQueue.addDataset(self.activeClipboard)
        return
    
    def _shutdownPipeline(self):
        # FIXME: this is a hack!!!!
        import lsst.ctrl.events as events
        import lsst.daf.base as dafBase

        shutdownTopic = "quitDayMOPS"      # hard-coded I think, as Stage does not have Pipeline policy
        externalEventTransmitter = events.EventTransmitter(self._evbroker, 
                                                           shutdownTopic)
        root = dafBase.PropertySet()
        externalEventTransmitter.publish(root)
        return



















"""
Simulate the flow of data into DiaSourceIDTonight

Used for testing purposes only.
"""
from DayMOPSStage import DayMOPSStage
import lib

import lsst.daf.persistence as persistence
from SafeDbStorage import SafeDbStorage





class TelescopeSimulatorStage(DayMOPSStage):
    """
    Simulate the telescope observing for one full night and populate the 
    DiaSourceIDTonight database table with the appropriate DiaSource IDs.
    
    This Stage should not be used in operations.
    """
    def __init__(self, stageId=-1, policy=None, startFromNightNumber=0):
        """
        Standard Stage initializer.
        """
        print "Entering TelescopeSimulatorStage.__init__!"
        super(TelescopeSimulatorStage, self).__init__(stageId, policy)
        
        # Get database details from policy
        self.dbLocStr = self.getValueFromPolicy('database')
        
        # Times in .paf files are always in days, we need hours.
        self.utOffset = self.getValueFromPolicy('utOffset', 
                                                lib.DEFAULT_UT_OFFSET) * 24.
        
        
        
        # Fetch the oldest DiaSource time (in UT TAI nsec).
        # Connect to the database.
        db = SafeDbStorage()
        db.setRetrieveLocation(persistence.LogicalLocation(self.dbLocStr))
        
        # If we do not have data between those two times, go to the next 
        # available times.
        # sql: select taiMidPoint from DiaSource order by taiMidPoint limit 1;
        db.setTableForQuery('DiaSource')
        db.outColumn('min(exposureStartTime)', True)      # isExpr=True
        db.outColumn('max(exposureStartTime)', True)      # isExpr=True
        db.query()
        if(db.next()):
            tMin = db.getColumnByPosDouble(0)
            tMax = db.getColumnByPosDouble(1)
            db.finishQuery()
        else:
            db.finishQuery()
            raise(Exception('No data in %s/DiaSource' %(self.dbLocStr, )))
        del(db)
        
        print "got tMin, tMax = ", tMin, tMax
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
        Copy all the DiaSources from night to night into the DiaSourceIDTonight
        table. Start from the oldest night not grater than 
        self.lastProcessedNight.
        """
        print "Entering TelescopeSimulatorStage.preprocess!"
        # Call the superclass preprocess() method.
        super(TelescopeSimulatorStage, self).preprocess()
        
        # Update the night counter.
        self.lastProcessedNight += 1
        
        # Wipe the table clean.
        db = SafeDbStorage()
        db.setPersistLocation(persistence.LogicalLocation(self.dbLocStr))
        db.startTransaction()
        db.executeSql('delete from DiaSourceIDTonight')
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
        
        # Copy the corresponding DiaSource IDs into DiaSourceIDTonight.
        # sql: insert into DiaSourceIDTonight (DiaSourceId) \
        #      (select diaSourceId from DiaSource \
        #       where exposureStartTime between %f and %f);
        sql = 'insert into DiaSourceIDTonight (DiaSourceId) \
               (select diaSourceId from DiaSource \
                where exposureStartTime between %f and %f)'
        db = SafeDbStorage()
        db.setRetrieveLocation(persistence.LogicalLocation(self.dbLocStr))
        db.startTransaction()
        db.executeSql(sql %(tMin, tMax))
        db.endTransaction()
        del(db)
        self.logIt('INFO', 'Inserted DiaSource IDs into DiaSourceIDTonight')
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



















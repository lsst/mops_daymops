"""
BASIC COURSE
System gets
1. a list of (NOT ATTRIBUTED AND NOT LINKED and NOT PRECOVERED) DIASources 
   belonging to last night.
as input.

A. System determines last night night number.

B. System puts the night number on the clipboard.


Policy
  1. Database name and location

Input
  1. [s for s in DIASources if (not s.ATTRIBUTED and not s.PRECOVERED and not 
      s.LINKED]

Output
  1. None/error code?

Side Effects
  DB Inserts
    1. None
  DB Updates
    1. None
  DB Deletes
    1. None
"""
from DayMOPSStage import DayMOPSStage
import DiaSourceList
import lib

import time




class SetupStage(DayMOPSStage):
    """
    This just computes the night number and puts it on the clipboard. Since 
    there is no communication between the non-parallel and the parallel part of 
    a Stage, we duplicate the computation in preprocess() and process().
    """
    def __init__(self, stageId=-1, policy=None):
        """
        Standard Stage initializer.
        """
        super(SetupStage, self).__init__(stageId, policy)
        
        # Read the configuration from policy.
        self.dbLocStr = self.getValueFromPolicy('database')
        self.utOffset = self.getValueFromPolicy('utOffset')
        return
    
    def preprocess(self):
        """
        Execute the non-parallel processing for the Intra-night linking Stage.
        
        Pseudo-code:
        1. Fetch last night's DIASources (at most one).
        2. Compute night number
        3. Update clipboard.
        """
        # Call the superclass preprocess.
        super(SetupStage, self).preprocess()
        self.logIt('INFO', 'Started preprocessing.')
        
        # Retrieve the DIASources to process.
        iter = DiaSourceList.diaSourceListForTonight(self.dbLocStr)
        try:
            source = iter.next()
        except StopIteration:
            source = None
        if(source):
            self.logIt('INFO', 'Found at least one DiaSource.')
        else:
            self.activeClipboard.put('nightNumber', None)
            self.logIt('INFO', 'Found 0 DiaSources.')
            return
        
        # Compute the night number.
        nightNumber = lib.mjdToNightNumber(source.getTaiMidPoint(),
                                           utOffset=self.utOffset)
        
        # Update the clipboard.
        self.activeClipboard.put('nightNumber', nightNumber)
        self.logIt('INFO', 'Updated the clipboard.')
        return
    
    def process(self):
        """
        Execute the non-parallel processing for the Intra-night linking Stage.
        
        Pseudo-code:
        1. Fetch last night's DIASources (at most one).
        2. Compute night number
        3. Update clipboard.
        """
        # Fetch the clipboard.
        self.activeClipboard = self.inputQueue.getNextDataset()
        self.logIt('INFO', 'Started processing.')
        
        # Retrieve the DIASources to process.
        iter = DiaSourceList.diaSourceListForTonight(self.dbLocStr)
        try:
            source = iter.next()
        except StopIteration:
            source = None
        if(source):
            self.logIt('INFO', 'Found at least one DiaSource.')
        else:
            self.activeClipboard.put('nightNumber', None)
            self.outputQueue.addDataset(self.activeClipboard)
            self.logIt('INFO', 'Found 0 DiaSources.')
            return
        
        # Compute the night number.
        nightNumber = lib.mjdToNightNumber(source.getTaiMidPoint(),
                                           utOffset=self.utOffset)
        
        # Update the clipboard.
        self.activeClipboard.put('nightNumber', nightNumber)
        
        # Put the clipboard back.
        self.outputQueue.addDataset(self.activeClipboard)
        self.logIt('INFO', 'Updated the clipboard.')
        return


    

    
    


















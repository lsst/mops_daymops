"""
Not implemented yet: just a basic skeleton


BASIC SCENARIO

System gets
1. a list of (NEW OR ATTRIBUTED) MovingObjects
2. a list of PROCESSED ScienceCCDExposures (i.e. RA, Dec, Filter, MJD, 
   FoVRadius)

A. System invokes "Compute Orbit-Field Intersection" on inputs 1. and 2.

B. For each MovingObject returned by step A.:
  B.1. System invokes "Compute Precise Ephemerides" on the MOvingObject and the 
       MJD of the Visits returned for that MovingObject by step A.

C. System groups MJDs in Input 2. by night.

D. For each night in step C.:
  D.1. System computes the union of the (RA, Dec, mag, MJD, ErrorEllipse) 
       returned by step B.1.
  D.2. System fetches the list of (NOT ATTRIBUTED AND NOT LINKED) Tracklets of 
       each Visit in input 2.
  D.3. System invokes "Identify Best DiaSource Match" on Ephemeris list from 
       step D.1. and list of Tracklets from step D.2.
  D.4. System invokes "Orbit Determination" on MovingObjects and candidate 
       Tracklets from step D.3.
  D.5. System stores the modified Orbits in the MovingObjects from step D.4.
  D.6. System markes the MovingObjects from step D.5. as NEW.

E. System marks MovigObjects from input 1. as PRECOVERED.

F. System markes appropriate DiaSources as PRECOVERED.



Policy
  1. MPC observatory OBSCODE
  2. FieldProximity config
  3. OrbitDetermination config
  4. Database name and location

Input
  1. a list of (NEW OR ATTRIBUTED) MovingObjects
  2. a list of PROCESSED ScienceCCDExposures (i.e. RA, Dec, Filter, MJD, 
     FoVRadius)
  3. the list of (NOT ATTRIBUTED AND NOT LINKED) Tracklets of each element in 
     input 2.

Output
  1. None/error code?

Side Effects
  DB Inserts
    1. None
  DB Updates
    1. Tracklet status -> LINKED
    2. MovingOject status -> PRECOVERED
  DB Deletes
    1. None
"""
from DayMOPSStage import DayMOPSStage
import attribution
import lib

import time




class PrecoveryStage(DayMOPSStage):
    def __init__(self, stageId=-1, policy=None):
        """
        Standard Stage initializer.
        """
        super(PrecoveryStage, self).__init__(stageId, policy)
        
        # Read the configuration from policy.
        self.obsCode = self.getValueFromPolicy('obsCode')
        self.residThreshold = self.getValueFromPolicy('residThreshold')
        self.maxSearchRadius = self.getValueFromPolicy('maxSearchRadius')
        self.minSearchRadius = self.getValueFromPolicy('minSearchRadius')
        self.maxArclengthForIod = self.getValueFromPolicy('maxArclengthForIod')
        self.uncertaintySigma = self.getValueFromPolicy('uncertaintySigma')
        self.dbLocStr = self.getValueFromPolicy('database')
        return
    
    def preprocess(self):
        super(PrecoveryStage, self).preprocess()
        self.outputQueue.addDataset(self.activeClipboard)
        return
    
    def process(self):
        # Fetch the clipboard.
        self.activeClipboard = self.inputQueue.getNextDataset()
        
        # What is our rank?
        i = self.getRank()
        n = self.getUniverseSize() - 1  # want only real slices
        self.logIt('INFO', 'Slice ID: %d/%d' %(i, n))
        
        
        # Do something!
        
        
        # Put the clipboard back.
        self.outputQueue.addDataset(self.activeClipboard)
        return
    
    
    


    

    
    


















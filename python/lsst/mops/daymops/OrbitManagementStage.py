"""
BASIC SCENARIO (not fully implemented yet).

System gets
1. a list of NEW MovingObjects
2. the list of (NOT NEW) MovingObjects
as input.

A. System invokes "Orbit Proximity"

B. For each MovingObject pair identified above:
  B.1. System retrieves the list of DIASources associated to each MovingObject.
  B.2. System invokes "Orbit Determination" on union of the two lists of 
       DIAsources.
  B.3. System marks the two MovingObjects as MERGED and creates a new 
       MovingObject linked to the 2 MERGED ones (using the new orbit derived 
       from the unified list of DIASources) marked as NEW.

RAINY DAY SCENARIO
If step B.2. (OrbitDetermination) failes, the two corresponding MovingObjects 
are left as they are (i.e. their status flags are not modified) and the loop 
continues.



Policy
  1. MPC observatory OBSCODE
  2. OrbitProximity config
  3. OrbitDetermination config
  4. Database name and location

Input
  1. a list of (NEW OR ATTRIBUTED OR PRECOVERED) MovingObjects
  2. a list of NOT (NEW OR ATTRIBUTED OR PRECOVERED) MovingObjects

Output
  1. None/error code?
"""
from DayMOPSStage import DayMOPSStage
from MovingObjectList import getAllPreliminaryMovingObjects, updateStatus
from MovingObject import STATUS as MO_STATUS
from Tracklet import STATUS as T_STATUS
import linking




class OrbitManagementStage(DayMOPSStage):
    def __init__(self, stageId=-1, policy=None):
        """
        Standard Stage initializer.
        """
        super(OrbitManagementStage, self).__init__(stageId, policy)
        
        # Read the configuration from policy.
        self.dbLocStr = self.getValueFromPolicy('database')
        return
    
    def preprocess(self):
        """
        Fetch all preliminary MovingObject instances from the database. Then, if
        >1 of these share the same Tracklet, choose the "best" one.
        
        Of course the meaning f best here is controversial (at the time of 
        writing, at least). The decision is made in linking.chooseOrbit().
        """
        # Call the superclass preprocess.
        super(OrbitManagementStage, self).preprocess()
        self.logIt('INFO', 'Started preprocessing.')
        
        # Get the prelim. MovingObject instances.
        movingObjectsIter = getAllPreliminaryMovingObjects(self.dbLocStr, 
                                                           shallow=False)
        
        # Now, get the full list of Tracklets and start examining each one to 
        # see if it is associated to >1 MovingObject. Not the best way to go 
        # about it since this would be a single SQL statement...
        # TODO: Implement something better here!
        trackletToMO = {}       # {tracklet: [MovingObject1,MovingObject2,...]}
        for mo in movingObjectsIter:
            tracklets = mo.getTracklets()
            for t in tracklets:
                tId = t.getTrackletId()
                if(trackletToMO.has_key(tId)):
                    trackletToMO[tId].append(mo)
                else:
                    trackletToMO[tId] = [mo, ]
        
        # And now back again!
        movingObjects = []
        for tId in trackletToMO.keys():
            # Do we have >1 MovingObject?
            if(len(trackletToMO[tId]) <= 1):
                # Nope!
                continue
            
            # Yup!
            mo = linking.chooseOrbit(trackletToMO[tId])
            if(not mo in movingObjects):
                # Update the status of this MovingObject and its Tracklets.
                mo.setStatus(MO_STATUS['NEW'])
                [t.setStatus(T_STATUS['ATTRIBUTED']) for t in mo.getTracklets()]
                movingObjects.append(mo)
        
        # Now update the status of the MovingObjects we have chosen and their
        # Tracklets as well.
        updateStatus(self.dbLocStr, movingObjects, updateTrackletStatus=True)
        return
    
    
    


    

    
    


















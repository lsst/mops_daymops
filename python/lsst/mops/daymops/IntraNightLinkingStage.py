"""
BASIC COURSE
System gets
1. a list of (NOT ATTRIBUTED AND NOT LINKED and NOT PRECOVERED) DIASources 
   belonging to the same night.
as input.

A. System invokes "Find tracklets"

B. System returns the newly formed linkages (a.k.a. Tracklets).

C. System flags the DIASources linked in each Tracklet as LINKED.

D. System flags the corresponding Visits as PROCESSED.


Policy
  1. FindTracklets config
  2. Database name and location

Input
  1. [s for s in DIASources if (not s.ATTRIBUTED and not s.PRECOVERED and not 
      s.LINKED]

Output
  1. None/error code?

Side Effects
  DB Inserts
    1. New Tracklets
  DB Updates
    1. DIASource statusm -> LINKED
  DB Deletes
    1. None
"""
from DayMOPSStage import DayMOPSStage
import DiaSourceList
import TrackletList
import linking




class IntraNightLinkingStage(DayMOPSStage):
    def __init__(self, stageId=-1, policy=None):
        """
        Standard Stage initializer.
        """
        super(IntraNightLinkingStage, self).__init__(stageId, policy)
        
        # Read the configuration from policy.
        self.maxV = self.getValueFromPolicy('maxV', linking.DEFAULT_MAXV)
        self.minObs = self.getValueFromPolicy('minObs', linking.DEFAULT_MINOBS)
        self.extended = self.getValueFromPolicy('extended', 
                                                linking.DEFAULT_EXTENDED)
        self.maxT = self.getValueFromPolicy('maxT', linking.DEFAULT_MAXT)
        self.collapseArgs = self.getValueFromPolicy('collapseArgs', 
                                                  linking.DEFAULT_COLLAPSE_ARGS)
        self.dbLocStr = self.getValueFromPolicy('database')
        return
    
    def preprocess(self):
        """
        Execute the non-parallel processing for the Intra-night linking Stage.
        
        Pseudo-code:
        1. Fetch last night's DIASources.
        2. Execute auton.findTracklets and get tracklets out.
        3. Update DB set linked DIASources.mopsStatus='L'
        4. insert into DB mops_Tracklet
        5. insert new entries into mops_TrackletsToDIASource.
        """
        # Call the superclass preprocess.
        super(IntraNightLinkingStage, self).preprocess()
        
        self.logIt('INFO', 'Starting processing.')
        
        # Retrieve the DIASources to process.
        sources = []
        for s in DiaSourceList.diaSourceListForTonight(self.dbLocStr):
            sources.append(s)
        self.logIt('INFO', 'Fetched %d DiaSources.' %(len(sources)))
        
        # Build tracklets out of the DiaSources we just selected. What we get 
        # out is a list of Tracklets each one of which lists the diaSourceIds
        # used to make that Tracklet.
        if(not sources):
            return
        
        # Build the tracklets.
        tracklets = [t for t in linking.trackletsFromDiaSources(sources)]
        
        # Put those Tracklets in the database.
        if(tracklets and len(tracklets)):
            self.logIt('INFO', 'Found %d Tracklets from %d DIASources.' \
                       %(len(tracklets), len(sources)))
            TrackletList.save(self.dbLocStr, tracklets)
        else:
            self.logIt('INFO', 'Found 0 Tracklets from %d DIASources.' \
                       %(len(sources)))
        return


    

    
    


















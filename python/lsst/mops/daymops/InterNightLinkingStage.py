"""
BASIC COURSE
System gets
1. a list of (NOT ATTRIBUTED AND NOT PRECOVERED AND NOT LINKED) Tracklets 
   belonging to >=1 nights within the last 30 days
as input.

A. System invokes "Link Tracklets"

B. System stores clusters of trackletIds, one cluster per Track


Policy
  1. LinkTracklets config
  3. MPC observatory OBSCODE
  4. Database name and location

Input
  1. a list of (NOT ATTRIBUTED AND NOT PRECOVERED AND NOT LINKED) Tracklets 
     belonging to >=1 nights within the last 30 days

Output
  1. None/error code?
"""
from DayMOPSStage import DayMOPSStage
import TrackletList
import TrackList
from Track import Track
from Tracklet import Tracklet
import linking
import lib




class InterNightLinkingStage(DayMOPSStage):
    def __init__(self, stageId=-1, policy=None):
        """
        Standard Stage initializer.
        """
        super(InterNightLinkingStage, self).__init__(stageId, policy)
        
        # Read the configuration from policy.
        self.obsCode = str(self.getValueFromPolicy('obsCode'))
        self.dbLocStr = self.getValueFromPolicy('database')
        
        self.slowMinV = self.getValueFromPolicy('slowMinV')
        self.slowMaxV = self.getValueFromPolicy('slowMaxV')
        self.slowVtreeThresh = self.getValueFromPolicy('slowVtreeThresh')
        self.slowPredThresh = self.getValueFromPolicy('slowPredThresh')
        self.slowMaxAccRa = self.getValueFromPolicy('slowMaxAccRa')
        self.slowMaxAccDec = self.getValueFromPolicy('slowMaxAccDec')
        
        self.fastMinV = self.getValueFromPolicy('fastMinV')
        self.fastMaxV = self.getValueFromPolicy('fastMaxV')
        self.fastVtreeThresh = self.getValueFromPolicy('fastVtreeThresh')
        self.fastPredThresh = self.getValueFromPolicy('fastPredThresh')
        self.fastMaxAccRa = self.getValueFromPolicy('fastMaxAccRa')
        self.fastMaxAccDec = self.getValueFromPolicy('fastMaxAccDec')
        
        self.plateWidth = self.getValueFromPolicy('plateWidth')
        
        self.minNights = self.getValueFromPolicy('minNights')
        self.timeSpan = self.getValueFromPolicy('timeSpan')
        self.utOffset = self.getValueFromPolicy('utOffset')        

        return
    
    def preprocess(self):
        """
        Create Track instances from Tracklets which are not part of any orbit 
        and have at least one DiaSource observed now - self.timeSpan days ago.
        """
        # We create Tracks in a preprocess() method rather than in the parallel
        # process since
        # 1. We do not have sky tessellation abilities.
        # 2. We do nopt want to lose Tracks just because we arbitrarily 
        #    partition Tracklets.
        
        # Call the superclass preprocess.
        super(InterNightLinkingStage, self).preprocess()
        self.logIt('INFO', 'Started preprocessing.')
        
        # Get the night number from the clipboard.
        nightNumber = self.activeClipboard.get('nightNumber')
        if(nightNumber == None):
            # No nightNumber on the clipboard. This is either an error or a sign
            # that there were no DiaSources for tonight. Either way say that we
            # quit and quit.
            self.logIt('INFO', 'No nightNumber on the clipboard. Quitting.')
            self.outputQueue.addDataset(self.activeClipboard)
            return
        self.logIt('INFO', 'Processing up to nightNumber %d' %(nightNumber))
        
        # From nightNumber and self.timeSpan, derive the minimum MJD to 
        # consider.
        minMjd = lib.nightNumberToMjdRange(nightNumber, self.utOffset)[0] - \
                 self.timeSpan
        
        # Fetch the list of available Tracklets within the last 30 days (meaning
        # fetch all non linked/attributed etc. Tracklets with at least one
        # DiaSource whose taiMidPoint is within 30 days of tonight's MJD).
        tracklets = []
        for tracklet in TrackletList.newTracklets(self.dbLocStr, 
                                                  shallow=False,
                                                  fromMjd=minMjd,
                                                  toMjd=None):
            tracklets.append(tracklet)
        self.logIt('INFO', 'Found %d tracklets.' %(len(tracklets)))
        if(not tracklets):
            self.outputQueue.addDataset(self.activeClipboard)
            return
        
        # Create Tracks from those Tracklets. We do this (internally) in two 
        # steps: a first step for slow movers (as defined in the policy file) 
        # and one for fast movers (again as defined in the policy file).
        tracks = []
        i = 0
        self.logIt('INFO', 'calling linkTracklets with slowMaxAccRa = %f, slowMaxAccDec = %f.'\
                       % (self.slowMaxAccRa, self.slowMaxAccDec))
        for rawTrack in linking.linkTracklets(tracklets, 
                                              self.slowMinV,
                                              self.slowMaxV,
                                              self.slowVtreeThresh,
                                              self.slowPredThresh,
                                              self.slowMaxAccRa,
                                              self.slowMaxAccDec,
                                              self.fastMinV,
                                              self.fastMaxV,
                                              self.fastVtreeThresh,
                                              self.fastPredThresh,
                                              self.fastMaxAccRa,
                                              self.fastMaxAccDec,
                                              self.minNights,
                                              self.plateWidth,
                                              self.obsCode):
            # We do not really need the full Tracklet object: just an instance
            # with the appropriate trackletId.
            i += 1
            tklts = [Tracklet(trackletId=tId) for tId in rawTrack]
            tracks.append(Track(trackId=None, tracklets=tklts))

        self.logIt('INFO', 'Found %d tracks.' %(len(tracks)))
        if(not tracks):
            self.outputQueue.addDataset(self.activeClipboard)
            return
        
        # Now just store the Tracks to the database but first wipe what is 
        # already there (since this is not parallel we are cool).
        TrackList.deleteAllTracks(self.dbLocStr)
        TrackList.save(self.dbLocStr, tracks)
        return
    
    
    


    

    
    


















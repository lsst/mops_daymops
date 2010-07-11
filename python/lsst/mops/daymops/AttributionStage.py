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
Not implemented yet: basic skeleton.


BASIC SCENARIO

System gets
1. a list of (NOT LINKED AND NOT ATTRIBUTED AND NOT PRECOVERED) Tracklets 
   belonging to the same night
2. the list of all NOT MERGED MovingObjects
3. the list of (RA, Dec, MJD, FieldRadius) of ScienceCCDExposures of the same 
   night as input 1.
4. the MPC OBSCODE of the observatory location
as input.

For each (RA, Dec, MJD, FieldRadius) in the input list 2:
  A. System invokes "Compute Orbit-Field of View Intersection" use case and 
     receives the appropriate subset of the input list of MovingObjects.
  B. For each Moving object in the list from step A:
    B.1. System invokes "Compute Precise Ephemeris" use case and receives the 
         list of predicted ephemeris (i.e. RA, Dec, mag, ErrorEllipse) for the 
         input MovingObject and MJD.
  C. System identifies most likely matches between the aggregated ephemeris list
     and the input list of Tracklets (by invoking "Identify Best DIASource 
     Match").
  D. System invokes "Orbit Determination" use case and recetives the list of 
     newly improved MovingObjects.
  E. System updates the system-level MovingObjects modified in step D (e.g. in 
     the database) with their new parameters.
  F. System markes MovingObjects from step E. as ATTRIBUTED.
  G. System flags the subset of the input Tracklets that were associated to 
     input MovingObjects in step F. as "LINKED".


Policy
  1. MPC observatory OBSCODE
  2. FieldProximity config
  3. OrbitDetermination config
  4. Database name and location

Input
  1. a list of (NOT LINKED AND NOT ATTRIBUTED AND NOT PRECOVERED) Tracklets 
     belonging to the same night
  2. the list of all NOT MERGED MovingObjects
  3. the list of (RA, Dec, MJD, FieldRadius) of ScienceCCDExposures of the same 
     night as input 1.

Output
  1. None/error code?

Side Effects
  DB Inserts
    1. None
  DB Updates
    1. Tracklet status -> LINKED
    2. MovingOject status -> ATTRINUTED
  DB Deletes
    1. None
"""
from DayMOPSStage import DayMOPSStage
from DiaSource import DiaSource
import TrackletList
import MovingObjectList
import attribution
import lib

import time




class AttributionStage(DayMOPSStage):
    def __init__(self, stageId=-1, policy=None):
        """
        Standard Stage initializer.
        """
        super(AttributionStage, self).__init__(stageId, policy)
        
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
        super(AttributionStage, self).preprocess()
        self.outputQueue.addDataset(self.activeClipboard)
        return
    
    def process(self):
        # Fetch the clipboard.
        self.activeClipboard = self.inputQueue.getNextDataset()
        
        # What is our rank?
        i = self.getRank()
        n = self.getUniverseSize() - 1  # want only real slices
        self.logIt('INFO', 'Slice ID: %d/%d' %(i, n))
        
        # Fetch all the  non attribbuted Tracklets with at least one DiaSource 
        # from tonight. DO not bother fetching the DiaSources for each Tracklet.
        tracklets = []
        for t in TrackletList.newTrackletsFromTonight(self.dbLocStr, 
                                                      shallow=True):
            tracklets.append(t)
        self.logIt('INFO', 'Found %d new Tracklets' %(len(tracklets)))
        if(not tracklets):
            self.logIt('INFO', 'Nothing to do for tonight.')
            self.outputQueue.addDataset(self.activeClipboard)
            return
        
        # We are going to come back to this Stage later since we need to create
        # orbits first.
        if(True):
            self.logIt('INFO', 'We\'ll come back to this Stage later.')
            self.outputQueue.addDataset(self.activeClipboard)
            return
            
        # Now, fetch 1/n of the known MovingObjects. We can do this because the 
        # search if on a per MovingObject basis.
        numObj = 0
        fetcher = MovingObjectList.getAllUnstableMovingObjects
        for obj in fetcher(self.dbLocStr, shallow=True, sliceId=i, numSlices=n):
            numObj += 1
            
            # Look for a tracklet that could be associated to obj.
            candidates = attribution.findCandidateTracklets(obj, tracklets)
            
            # 
            
            
        self.logIt('INFO', 'Found %d MovingObjects' %(numObj))
        
        # Put the clipboard back.
        self.outputQueue.addDataset(self.activeClipboard)
        return
    
    
    


    

    
    


















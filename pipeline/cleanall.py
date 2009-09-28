#!/usr/bin/env python
import sys
import lsst.daf.persistence as persistence


if(len(sys.argv) != 2):
    sys.stderr.write('usage: cleanall.py <database name>\n')
    sys.stderr.flush()
    sys.exit(1)
database = sys.argv[1]

db = persistence.DbStorage()
db.setRetrieveLocation(persistence.LogicalLocation('mysql://localhost:3306/%s' %(database)))
db.startTransaction()
db.executeSql('delete from DIASourceIDTonight')
db.executeSql('delete from mops_Tracklet')
db.executeSql('delete from mops_TrackletsToDIASource')
db.executeSql('delete from mops_TracksToTracklet')
db.endTransaction()
sys.exit(0)

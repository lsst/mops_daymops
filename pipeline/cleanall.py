#!/usr/bin/env python

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

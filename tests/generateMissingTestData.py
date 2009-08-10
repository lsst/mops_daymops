import random
import MySQLdb as DBI


dbh = DBI.connect(db='mops_test', user='fpierfed', passwd='f1gata')
dbc = dbh.cursor()




if(False):
    sql = 'select movingObjectId from MovingObject'
    n = dbc.execute(sql)
    rows = dbc.fetchall()
    ids = zip(*rows)[0]
    n = len(ids)
    
    # We want 1000 moving objects!
    chosenIds = []
    for i in range(1000):
        _id = str(int(random.random() * n))
        while(_id in chosenIds):
            _id = str(int(random.random() * n))
        chosenIds.append(_id)
    
    # remove everything else.
    sql = 'delete from MovingObject where movingObjectId not in (%s)'
    n = dbc.execute(sql %(','.join(chosenIds)))
if(False):
    # generate fake MovingObject -> Tracklet info.
    sql = 'select movingObjectId from MovingObject'
    n = dbc.execute(sql)
    rows = dbc.fetchall()
    moIds = zip(*rows)[0]
    noTot = len(moIds)
    
    sql = 'select trackletId from mops_Tracklet'
    n = dbc.execute(sql)
    rows = dbc.fetchall()
    tIds = zip(*rows)[0]
    tTot = len(tIds)
    
    # we want each MovingObject to have between 0 and 10 tracklets. Overlaps are
    # fine.
    sql = 'insert into mops_MovingObjectToTracklet (movingObjectId, trackletId)'
    sql += ' values '
    for moId in moIds:
        valSQL = ''
        n = int(random.uniform(0, 10))
        if(not n):
            continue
        # choose n random Tracklets.
        idxs = []
        for i in range(n):
            idx = int(random.random() * (tTot - 1))
            while(idx in idxs):
                idx = int(random.random() * (tTot - 1))
            idxs.append(idx)
            valSQL += '(%d, %d), ' %(moId, tIds[idx])
        valSQL = valSQL[:-2]
        
        # Insert the values.
        # print(sql + valSQL)
        dbc.execute(sql + valSQL)
if(False):
    # Generate random Track linkages.
    sql = 'select trackletId from mops_Tracklet'
    n = dbc.execute(sql)
    rows = dbc.fetchall()
    trIds = zip(*rows)[0]
    trTot = len(trIds)
    
    # Generate 1000 trackIds
    tIds = range(1, 1001, 1)
    tTot = len(tIds)
    
    # we want each Track to have between 0 and 10 tracklets. Overlaps are fine.
    sql = 'insert into mops_TracksToTracklet (trackId, trackletId) values '
    for tId in tIds:
        valSQL = ''
        n = int(random.uniform(0, 10))
        if(not n):
            continue
        # choose n random Tracklets.
        idxs = []
        for i in range(n):
            idx = int(random.uniform(1, trTot))
            while(idx in idxs):
                idx = int(random.uniform(1, trTot))
            idxs.append(idx)
            valSQL += '(%d, %d), ' %(tId, trIds[idx])
        valSQL = valSQL[:-2]
        
        # Insert the values.
        # print(sql + valSQL)
        dbc.execute(sql + valSQL)























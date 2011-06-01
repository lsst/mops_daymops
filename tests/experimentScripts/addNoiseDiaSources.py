#!/usr/bin/env python

""" jmyers Oct 5 2010

Tools for adding noise detections randomly throughout simulated image catalogs.

Note on image schedules, since I always forget this myself:

opsim holds the schedule of images observed.  Each image is uniquely
identified by its obsHistId (I think), and this is the data which
diaSources hold.  Each image and falls in a given field. All visits to
a given field should have the same center RA and Dec CURRENTLY but in
the long run they will not.

Your input/output diaSourceTable should look like this:

# +-------------+------------+------+-----+---------+-------+
# | Field       | Type       | Null | Key | Default | Extra |
# +-------------+------------+------+-----+---------+-------+
# | diaSourceId | bigint(20) | YES  | MUL | NULL    |       | 
# | opSimId     | bigint(20) | YES  |     | NULL    |       | 
# | ssmId       | bigint(20) | YES  |     | NULL    |       | 
# | ra          | double     | YES  |     | NULL    |       | 
# | decl        | double     | YES  |     | NULL    |       | 
# | taiMidPoint | double     | YES  |     | NULL    |       | 
# | mag         | double     | YES  |     | NULL    |       | 
# | snr         | double     | YES  |     | NULL    |       | 
# +-------------+------------+------+-----+---------+-------+


"""




import MySQLdb as db
import numpy 
import rotateCoords 
import random
import sys

IMAGE_RADIUS_DEG=1.784


DB_HOST="localhost"
DB_USER="jmyers"
DB_PASS="jmyers"
OPSIM_DB="opsim_3_61"
OPSIM_TABLE="output_opsim3_61"

FALSE_DIA_SSMID=-1
FALSE_DIA_MAG=1
FALSE_DIA_SNR=1



def getAllFieldsObserved(cursor, trueDiasDb, trueDiasTable):
    """ Return list of all obsHistIds for all fields visited in our
    simulation. """
    # TBD: Currently we just check for fields which contained some
    # diaSource; this could be a problem if we visited a field but it
    # had no diasources.  This is good enough for now since we are
    # looking at the ecliptic which is very dense but in the general
    # case we may need to do something more intelligent.
    

    sql = """select distinct(opSimId) from %s.%s;""" % (trueDiasDb, trueDiasTable)
    cursor.execute(sql)
    obsHists = cursor.fetchall()
    # we get a list of singleton lists: [[1],[2],[3]...]. Change it.
    obsHists = map (lambda x: x[0], obsHists)
    return obsHists



def getImageTimeLoc(cursor, obsHistId):
    """ return (image center ra, dec, image Time) for a given obsHistId. Return RA, Dec in DEGREES. """
    sql = """select fieldRA, fieldDec, expMJD from %s.%s where obsHistId=%d;""" % (OPSIM_DB, OPSIM_TABLE, obsHistId)
    cursor.execute(sql)
    [ra,dec,time] = cursor.fetchone()
    return numpy.degrees(ra), numpy.degrees(dec), time




def addFalseDias(locs, time, obsHistId, firstId, cursor,
                 falseDiasDb, falseDiasTable):
    """ Add false dias with specified locations. They are assumed to
    share a common time and obsHistId.  Ids are assigned starting at
    firstId and climbing upward. All other needed values are taken
    from the constants defined in this file"""

    s = """ INSERT INTO %s.%s (diaSourceId, opSimId, ssmId, 
              ra, decl, taiMidPoint, mag, snr) VALUES """ % \
        (falseDiasDb, falseDiasTable)
    nextId = firstId
    for i in range(len(locs)):
        loc = locs[i]
        s += "(%d, %d, %d, %12f, %12f, %12f, %12f, %12f)" % \
            (nextId, obsHistId, FALSE_DIA_SSMID, loc[0], loc[1],
             time, FALSE_DIA_MAG, FALSE_DIA_SNR)
        if i != len(locs) - 1:
            s += ", "
        nextId += 1
    s += ";"
    #print s
    cursor.execute(s)



def angularDistance(a, b):
    """ return distance between a and b, where a and b are angles in degrees. """
    while abs(a - b) > 180:
        if a > b:
            b += 360.
        else:
            a += 360.
    return a - b

def convertToStandardDegrees(angle):
    while angle > 360.:
        angle -= 360.
    while angle < 0.:
        angle += 360.
    return angle


def greatCircleDistance(RA0, Dec0, RA1, Dec1):
    """
    return the great-circle distance between two points on the sky,
    uses haversine formula
    """
    deg_to_rad = numpy.pi / 180.
    rad_to_deg = 180. / numpy.pi

    RADist = angularDistance(RA0, RA1);
    DecDist = angularDistance(Dec0, Dec1);    
    #convert all factors to radians
    RADist = deg_to_rad*convertToStandardDegrees(RADist);
    DecDist = deg_to_rad*convertToStandardDegrees(DecDist);
    Dec0 = deg_to_rad*convertToStandardDegrees(Dec0);
    Dec1 = deg_to_rad*convertToStandardDegrees(Dec1);
    r = 2*numpy.arcsin(numpy.sqrt( (numpy.sin(DecDist/2.))**2 + 
                                 numpy.cos(Dec0)*numpy.cos(Dec1)*(numpy.sin(RADist/2))**2) );
    #back to degrees
    return rad_to_deg*r;





def chooseRandomLocationsInImage(n, centerRa, centerDec, radius):
    """
    Chooses n random locations within radius of centerRa, centerDec
    (all in DEGREES) and returns their RA, Dec (in DEGREES)
    """
    chosenPoints = 0
    pointsRa = []
    pointsDec = []
    #choose points < radius away from the origin.
    while chosenPoints < n:
        ra = random.random() * (radius * 2) - radius
        dec = random.random() * (radius * 2) - radius
        dist = greatCircleDistance(ra, dec, 0., 0.)
        if dist < radius:
            pointsRa.append(ra)
            pointsDec.append(dec)
            chosenPoints += 1
    
    #rotate coordinates from origin to center of image.
    newRas, newDecs =  rotateCoords.rotateMany(0, centerRa, 
                                               0, centerDec, 
                                               pointsRa, pointsDec)
    chosenPoints = []
    for i in range(len(newRas)):
        chosenPoints.append([newRas[i], newDecs[i]])

    # sanity check
    for i in range(len(chosenPoints)):
        pt = chosenPoints[i]
        [ra, dec] = pt
        dist = greatCircleDistance(ra, dec, centerRa, centerDec)
        if dist > radius:
            raise Exception("Got point %f %f was %f from center location\n" % (ra, dec, dist))
    return chosenPoints
    


def getLastTrueDiaId(cursor, trueDiasDb, trueDiasTable):
    """ return the maximum diaSource id in the diaSource table. """
    s = """ SELECT MAX(diaSourceId) FROM %s.%s; """ % \
        (trueDiasDb, trueDiasTable)
    cursor.execute(s)
    return cursor.fetchone()[0]


if __name__=="__main__":
    
    noisePointsPerImage = int(sys.argv[1])

    trueDiasDb = sys.argv[2]
    trueDiasTable = sys.argv[3]

    falseDiasDb = sys.argv[4]
    falseDiasTable = sys.argv[5]

    print "Copying dias from ", trueDiasDb , "." , trueDiasTable
    print "Adding dias and ", noisePointsPerImage, " noise points per image to table ", \
        falseDiasDb, ".", falseDiasTable

    conn = db.connect(host=DB_HOST, user=DB_USER, passwd=DB_PASS)
    curs = conn.cursor()

    #curs.execute("DROP TABLE IF EXISTS %s.%s" % (falseDiasDb, falseDiasTable))
    curs.execute("CREATE TABLE %s.%s LIKE %s.%s;" % \
                     (falseDiasDb, falseDiasTable,
                      trueDiasDb, trueDiasTable))

    curs.execute("INSERT INTO %s.%s SELECT * FROM  %s.%s;" % \
                     (falseDiasDb, falseDiasTable, 
                      trueDiasDb, trueDiasTable))
    
    obsHists = getAllFieldsObserved(curs, falseDiasDb, falseDiasTable)
    firstFalseDiaId = getLastTrueDiaId(curs, falseDiasDb, falseDiasTable) +  1
    curFalseDiaId = firstFalseDiaId
    for obsHist in obsHists:
        centerRa,centerDec,time = getImageTimeLoc(curs, obsHist)
        print centerRa, centerDec, time
        
        locs = chooseRandomLocationsInImage(noisePointsPerImage, 
                                            centerRa, centerDec, 
                                            IMAGE_RADIUS_DEG)

        addFalseDias(locs, time, obsHist,  curFalseDiaId, curs, 
                     falseDiasDb, falseDiasTable)
        curFalseDiaId += len(locs)
     

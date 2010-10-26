#!/usr/bin/env python

""" jmyers sep 1

Some python utilities for rotating RA, Dec points using a rotation
matrix and converting to/from Cartesian coords.

"""



from numpy import *



def toCartesian(ra, dec):
    """ return x, y, z cartesian points given ra, dec pair in DEGREES"""
    x = cos(radians(ra)) * sin(radians(90 - dec))
    y = sin(radians(ra)) * sin(radians(90 - dec))
    z = cos(radians(90 - dec))
    return x,y,z

def toRaDec(x, y, z):
    """ given a cartesian point in x,y,z, return the corresponding point in degrees."""
    r = degrees(arctan2(y, x))
    d = degrees(arctan2(z, sqrt(x*x + y*y)))
    return r,d

def buildRotationMatrix(raOld, raNew, decOld, decNew):
    """ build a rotation matrix which will rotate points from old RA,
    Dec (in degrees) to new RA, Dec in degrees. Note that this matrix
    is to be applied to CARTESIAN points, though!"""

    raOld, raNew, decOld, decNew = map(radians, [raOld, raNew, decOld, decNew])

    temp1 = array([[ cos(raOld), sin(raOld), 0 ],
                   [-sin(raOld), cos(raOld), 0 ],
                   [0,        0,       1 ]])

    dDec = decNew - decOld

    temp2 = array([[ cos(dDec),  0,  -sin(dDec) ],
                   [ 0,          1,  0          ],
                   [ sin(dDec),  0,  cos(dDec)  ]])

    temp3 = array([[ cos(-raNew), sin(-raNew), 0 ],
                   [-sin(-raNew), cos(-raNew), 0 ],
                   [0,        0,       1 ]])
    return dot(dot(temp3, temp2), temp1)


def rotateMany(raOld, raNew, decOld, decNew, ras, decs):
    """All units in degrees. given a set of RAs and Decs (as lists), rotate
    all detections from raOld to raNew and decOld to decNew """

    rot = buildRotationMatrix(raOld, raNew, decOld, decNew)
    newRas = []
    newDecs = []
    for i in range(len(ras)):
        ra, dec = ras[i], decs[i]
        x,y,z = toCartesian(ra, dec)
        [[xp], [yp], [zp]] = dot(rot, array([[x],
                                             [y],
                                             [z]]))    
        

        rap, decp = toRaDec(xp, yp, zp)
        newRas.append(rap)
        newDecs.append(decp)

    return newRas, newDecs



EPSILON=1e-10
def eq(a, b): 
    return abs(a - b) < EPSILON

import sys

if __name__ == "__main__":
    """ a few simple tests """
    ra, dec = 100, 30
    print ra, dec
    x,y,z = toCartesian(ra, dec)
    print x, y, z
    ra2, dec2 = toRaDec(x,y,z)
    print ra2, dec2
    if not (eq(ra,ra2) and eq(dec,dec2)):
        print "TEST FAILED"
        sys.exit(1)

    rot = buildRotationMatrix(ra, ra - 30, dec, dec - 20)
    [[xp], [yp], [zp]] = dot(rot, array([[x],
                                         [y],
                                         [z]]))    
    rap, decp = toRaDec(xp, yp, zp)
    print rap, decp
    if not (eq(rap, ra - 30) and eq(decp, dec - 20)):
        print "TEST FAILED"
        sys.exit(1)

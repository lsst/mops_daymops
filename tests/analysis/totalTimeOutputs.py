#!/usr/bin/env python

""" jmyers oct 25 2010

This is a script for taking a series of outputs from GNU time and summing up the total time, reporting it in seconds.

We expect your time to write output in the following format: 

---- start example
12261.42user 2.83system 3:24:27elapsed 99%CPU (0avgtext+0avgdata 0maxresident)k
309992inputs+6424outputs (6major+379772minor)pagefaults 0swaps
----- end example

"""



import glob
import sys


def getTimeInSecsFromFile(infile):
    lines = file(infile,'r').readlines()
    tokens = lines[0].split()
    if tokens[0][-4:] != "user" or \
            tokens[1][-6:] != "system" or\
            tokens[2][-7:] != "elapsed":
        print infile
        print tokens[0][-4:]
        print tokens[1][-6:]
        print tokens[2][-7:]
        raise Exception("Incorrectly formatted time file!?")
    
    user = float(tokens[0][:-4])
    sys = float(tokens[1][:-6])
    elapsedString = tokens[2][:-7]
    print elapsedString
    elapsedParts = elapsedString.split(':')
    if len(elapsedParts) == 3:
        h,m,s = map(float, elapsedParts)
    elif len(elapsedParts) == 2:
        h = 0
        m,s = map(float, elapsedParts)
    elif len(elapsedParts) == 1:
        h = 0
        m = 0
        s = float(elapsedParts[0])
    
    elapsed = h * 60 * 60 + m * 60 + s 
    print h, m, s, elapsed 

    return elapsed, user, sys


if __name__=="__main__":
    if len(sys.argv) != 2:
        print "USAGE: totalTimeOutputs timeOutputGlob"
        sys.exit(1)

    infiles = glob.glob(sys.argv[1])

    print "Reading files: ", infiles

    totalElapsed = 0
    totalUser = 0
    totalSys = 0

    for infile in infiles:
        elapsed, user, sys = getTimeInSecsFromFile(infile)
        totalElapsed += elapsed
        totalUser += user
        totalSys += sys

    print "Given ", len(infiles), " files, total CPU usage in seconds was: "
    print "Elapsed: ", totalElapsed, "\tUser: ", totalUser, "\tSys: ", totalSys

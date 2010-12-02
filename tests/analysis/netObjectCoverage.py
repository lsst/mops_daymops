#!/usr/bin/env python

""" jmyers dec 1 2010

Parse through the output files from evalObjectCoverage and evaluate
the net object coverage from many runs, weighted in a sane manner.

usage: netObjectCoverage.py <glob of all files>
"""



import sys
import glob


def getCoverageFromLines(lines):
    lines = map(lambda x: x.split(), lines)
    cov = float(lines[0][2])
    num = int(lines[1][2])
    return [cov, num]


def netObjectCoverageManyFiles(infileNames):
    netObjects = 0
    netCoverage = 0.0

    for inf in infileNames:
        #print inf
        lines = file(inf,'r').readlines()
        [coverage, numObj] = getCoverageFromLines(lines)
        netCoverage += coverage*numObj
        netObjects += numObj
    return [netCoverage / float(netObjects), netObjects]


if __name__=="__main__":
    infiles = glob.glob(sys.argv[1])
    print "looking in files: ", infiles

    [objectCoverage, numObjects] = netObjectCoverageManyFiles(infiles)
    print "Average object coverage: ", objectCoverage
    print "Number of objects: ", numObjects
        

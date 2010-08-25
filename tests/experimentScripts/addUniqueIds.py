#!/usr/bin/env python
""" 

jmyers may 11 2010

"""

import sys


if __name__=="__main__":
    inf = file(sys.argv[1], 'r')
    outf = file(sys.argv[2], 'w')
    
    line = inf.readline()
    curId = 0
    while line != "":
        items = line.split()
        outf.write(str(curId) + " " + " ".join(items) + '\n')
        curId += 1
        line = inf.readline()
    

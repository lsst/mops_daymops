#!/usr/bin/env python

""" jmyers

select a random subset of a series of items.  Set of items should fit
in memory.

Useful for getting e.g. 50%, 25% of the objects seen in our simulated
survey.

"""


import random


def selectSubset(items, nToChoose):
    """ randomly select and return nToChoose elements of the list list items."""
    outList = []
    indicesToKeep = set()
    nItems = len(items)
    # choose nItems unique indices
    for i in range(nToChoose):  
        choice = random.randint(0, nItems - 1)
        while choice in indicesToKeep:
            choice = random.randint(0, nItems - 1)
        indicesToKeep.add(choice)
    # add those items to outList.
    for i in indicesToKeep:
        outList.append(items[i])
    return outList


def writeItems(outfile, items):
    for i in items:
        outfile.write("%d\n" % i)
    

if __name__=="__main__":
    import sys
    items = map(int, file(sys.argv[1],'r').readlines())
    fraction = float(sys.argv[2])
    randomSeed = int(sys.argv[3])
    outfile = file(sys.argv[4],'w')
    print "Selecting ", fraction*100, "% of the items from ", sys.argv[1], \
        " and writing them to ", outfile
    print "Seeding random number generator with value ", randomSeed
    random.seed(randomSeed)

    nItems = len(items)

    selectedItems = selectSubset(items, int(nItems *fraction))

    writeItems(outfile, selectedItems)
    

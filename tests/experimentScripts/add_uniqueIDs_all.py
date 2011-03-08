#!/usr/bin/env python
""" 

ljones 
$Id$ 

Add unique diasource ID's to all sources in a series of diasource files.
(check format of file to see whether to edit this script)

# diasource ssmid obshistid filter expmjd ra ra_err dradt dec dec_err ddecdt
   dist mag mag_err snr 

"""

import sys
import os

if __name__=="__main__":
    """Give this script the root of the diasource filenames, it adds unique 
    diasourceID's to each file. """

    fileroot = sys.argv[1]
    tfilelist = os.listdir(".")
    infilelist = []
    for t in tfilelist:
        if t.startswith(fileroot):
            infilelist.append(t)
    
    diacount = 0 
    
    outfile = open("all_dias", "w")

    for infile in infilelist:
        f = open(infile, "r")
        for line in f:
            values = line.split()
            if values[0].startswith('diasourceid'):
                continue
            if values[0].startswith("#"):
                continue
            # assume values[0] is previous, non-unique diasourceID
            outfile.write("%s " %(diacount) + " ".join(values[1:]) + '\n')
            diacount = diacount + 1
        f.close()

# script to help generate diasources and their magnitudes /SNR in a particular filter

import os
import numpy
import selfcal.analyze.useful_input as ui
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass

# filter id table:
filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
filterdict = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'y':5}
rootdir = os.getenv("LSST_THROUGHPUTS_BASELINE")
lsst = {}
for f in filterlist:
    lsst[f] = Bandpass()
    lsst[f].readThroughput(rootdir + "/total_" + f + ".dat")
imsimband = Bandpass()
imsimband.imsimBandpass()

# set up colors of objects in all bandpasses
rootdir_sed = "."
sedlist = ('S.dat', 'C.dat')
simobjs = {}
simmags = {}
for sed in sedlist:
    simobjs[sed] = Sed()
    simobjs[sed].readSED_flambda(rootdir_sed + sed)
    simmags[sed] = {}
    for filter in filterlist:
        simmags[sed][filter] = simobjs[sed].calcMag(lsst[filter])
    simmags[sed]['ims'] = simobjs[sed].calcMag(imsimband)
    

# open and read file containing obshistids (key to catalog and to opsim)
infile = 'obsids'
f = open(infile, 'r')
outfile = 'obsinfo'

db,cursor = ui.sqlConnect()

printheader = True
newfile = True

diasourceid = 0

for line in f:
    (obshistid,) = line.split()
    sqlquery = "select obshistid, expmjd, filter, 5sigma_ps, night, fieldid from output_opsim3_61 where obshistid = %s and propid!=217" %(obshistid)
    results = ui.sqlQuery(cursor, sqlquery, verbose=False)
    keys = ('obshistid', 'expmjd', 'filter', '5sigma', 'night', 'fieldid')
    data = ui.assignResults(results, keys)
    ui.writeDatafile(outfile, data, keys, printheader=printheader, newfile=newfile)
    printheader = False
    newfile = False
    
    ampexposureid = obshistid

    catalogin = "catalog_" + obshistid + "_ssm_0"
    cin = open(catalogin, 'r')
    catalogout = "diasources_" + obshistid + "_%d" %(data['night'][0])
    print "Writing new catalog %s" %(catalogout)
    cout = open(catalogout, 'w')
    for detection in cin:
        # read value from catalog file
        values = detection.split()
        objid = values[1]
        ra = values[2]  #deg
        dec = values[3]  #deg
        imsimmag = float(values[4])
        sed = values[5][-5:]
        # calculate magnitude and SNR of this object in output
        dmag = imsimmag - simmags[sed]['ims'] 
        magout = simmags[sed][data['filter'][0]] + dmag
        flux_ratio = n.power(10, 0.4*(data['5sigma'][0]-magout))
        snr = 5 * (flux_ratio)
        if snr < 5: 
            #print "dumping object %s at %f mag, %f snr: background %f" %(objid, magout, snr, data['5sigma'])
            continue #dump this object if faint
        # calculate astrometric error in ra/dec
        rgamma = 0.039
        # average seeing is 0.7" (or 700 mas) 
        xval = n.power(10, 0.4*(magout - data['5sigma'][0]))
        error_rand = 0.7 * n.sqrt((0.04-rgamma)*xval + rgamma*xval)
        error_sys = 0.1 
        astrom_error = n.sqrt(error_sys * error_sys + error_rand*error_rand) / 60.0 / 60.0 # deg
        # calculate magnitude error         
        error_sys = 0.005     # systematic error - 0.005
        error_rand = n.sqrt((0.04-rgamma)*xval + rgamma*xval*xval)
        mag_error = n.sqrt(error_sys*error_sys + error_rand*error_rand)     
        ssmid = objid
        print >>cout, diasourceid, obshistid, ra, astrom_error, dec, astrom_error, expmjd, \
            mag, mag_error, SNR, filter,  
        # increment diasource id (unique across objects)
        diasourceid = diasourceid + 1
    # done with all observations    
    cin.close()
    cout.close()
# done with all obshistids
f.close()


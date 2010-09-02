import numpy
import useful_input as ui

def calcAstrometricError(mag, m5, seeing):
    """Calculate the astrometric error, for object catalog purposes.
    
    Returns astrometric error for a given SNR, in arcseconds"""
    rgamma = 0.039
    xval = numpy.power(10, 0.4*(mag-m5))
    # The average seeing is 0.7" (or 700 mas).
    # seeing is given in arcseconds, as in opsim
    error_rand = seeing* numpy.sqrt((0.04-rgamma)*xval + rgamma*xval*xval)
    # The systematic error floor in astrometry:  0.01 " ... but for now let's use a larger value
    # that is more conservative (and seems very believable) of 0.1"
    error_sys = 0.1
    astrom_error = numpy.sqrt(error_sys * error_sys + error_rand*error_rand)
    return astrom_error

if __name__ == '__main__':
    
    f = open('dias_pt1_nodeep.short', 'r')
    #f = open('test', 'r')
    fout = open('dias_pt1_nodeep.short.astromErr', 'w')
    prevopsimid = 0
    conn, cursor = ui.sqlConnect(hostname='lsst-db.astro.washington.edu', 
                                 username='lsst', passwdname='lsst', dbname='opsimdev')
    for line in f:
        values = line.split()
        opsimid = int(values[0])
        ssmid = int(values[1])
        ra = float(values[2])
        decl = float(values[3])
        expmjd = float(values[4])
        mag = float(values[5])
        snr = float(values[6])
        
        if opsimid != prevopsimid :
            query = "select 5sigma_ps, seeing from output_opsim3_61 where obsHistID=%d" %(opsimid)
            results = ui.sqlQuery(cursor, query, verbose=False)
            keys = ('mag5', 'seeing')
            data = ui.assignResults(results, keys)
            prevopsimid = opsimid
        
        astromErr = calcAstrometricError(mag, data['mag5'], data['seeing'])
        astromErr = astromErr / 3600.0
        
        ra_out = ra + numpy.random.normal(loc=0, scale=astromErr)
        if ra_out < 0 :
            ra_out = ra_out + 360.0
        if ra_out > 360 : 
            ra_out = ra_out % 360.0
        decl_out = decl + numpy.random.normal(loc=0, scale=astromErr)


        print >>fout, opsimid, ssmid, ra_out, decl_out, expmjd, mag, snr

    f.close()
    fout.close()
    ui.sqlEndConnect(conn, cursor)


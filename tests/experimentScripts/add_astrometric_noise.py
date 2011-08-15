import numpy
import useful_input as ui
import sys

def calcAstrometricError(mag, m5, seeing, systematic_error=0.1):
    """Calculate the astrometric error, for object catalog purposes.
    
    Returns astrometric error for a given SNR, in arcseconds"""
    rgamma = 0.039
    xval = numpy.power(10, 0.4*(mag-m5))
    # The average seeing is 0.7" (or 700 mas).
    # seeing is given in arcseconds, as in opsim
    error_rand = seeing* numpy.sqrt((0.04-rgamma)*xval + rgamma*xval*xval)
    # The systematic error floor in astrometry:  0.01 " ... but for now let's use a larger value
    # that is more conservative (and seems very believable) of 0.1"
    error_sys = systematic_error
    astrom_error = numpy.sqrt(error_sys * error_sys + error_rand*error_rand)
    return astrom_error

if __name__ == '__main__':
    
    if len(sys.argv)<2:
        print "usage: python add_astrometric_error.py INFILE OUTFILE"
        exit()
    infile = sys.argv[1]
    outfile = sys.argv[2]
    print "Reading from %s and writing to %s" %(infile, outfile)
    
    f = open(infile, 'r')
    #f = open('test', 'r')
    fout = open(outfile, 'w')
    conn, cursor = ui.sqlConnect(hostname='localhost', 
                                 username='jmyers', passwdname='jmyers', dbname='opsim_3_61')

    # ordered list of input data 
    inputdata_keys = ('diasourceID', 'opsimID', 'ssmID', 'ra', 'decl', 'expmjd', 'mag', 'snr')
    inputdata_keytypes = ('int', 'int', 'int', 'float', 'float', 'float', 'float', 'float')
    data = {} 
    prevopsimid = 0
    for line in f:
        values = line.split()
        for i in range(len(inputdata_keys)):
            if inputdata_keytypes[i] == 'int':
                data[inputdata_keys[i]] = int(values[i])
            else:
                data[inputdata_keys[i]] = float(values[i])

        if data['opsimID'] != prevopsimid :
            query = "select 5sigma_ps, seeing from output_opsim3_61 where obsHistID=%d" %(data['opsimID'])
            results = ui.sqlQuery(cursor, query, verbose=False)
            opsim_keys = ('mag5', 'seeing')
            opsim_data = ui.assignResults(results, opsim_keys)
            prevopsimid = data['opsimID']

        # calculate astrometric error in "
        astromErr = calcAstrometricError(data['mag'], opsim_data['mag5'], opsim_data['seeing'])
        # turn astrometric error into degrees to match ra/dec
        astromErr = astromErr / 3600.0
        
        data['ra'] = data['ra'] + numpy.random.normal(loc=0, scale=astromErr)
        if data['ra'] < 0 :
            data['ra'] = data['ra'] + 360.0
        if data['ra'] > 360 : 
            data['ra'] = data['ra'] % 360.0
        data['decl'] = data['decl'] + numpy.random.normal(loc=0, scale=astromErr)

        writestring = ""
        for i in range(len(inputdata_keys)):
            if inputdata_keytypes[i] == 'int':
                writestring = writestring + "%d " %(data[inputdata_keys[i]])
            else:
                if (inputdata_keys[i] == 'ra') | (inputdata_keys[i] == 'decl'):
                    writestring = writestring + "%.9f " %(data[inputdata_keys[i]])
                else: 
                    writestring = writestring + "%f " %(data[inputdata_keys[i]])
        print >>fout, writestring
        
    f.close()
    fout.close()
    ui.sqlEndConnect(conn, cursor)


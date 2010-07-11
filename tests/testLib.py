#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import math
import random
import unittest


try:
    import lsst.mops.daymops.lib as lib
except:
    raise(ImportError('Please setup daymops first.'))


# Ground truth.
ra = [11.701302, 98.015470, 232.439537, 143.395192, 24.881267, 293.523715,
      256.655854, 168.365205, 307.907172, 332.772905]
dec = [85.193105, -85.221559, -49.463942, -26.768707, -30.393225, -67.651531,
       7.824027, -36.874703, 67.380175, 11.579160]
# Disances beween, ra[0], dec[0] and ra[1], dec[1]; ra[0], dec[0] and ra[2], dec[2]
# all combination. Distances computed with code from
# http://adfwww.gsfc.nasa.gov/asca/docs/proc_script/rev1/doc/objdist.f
# Modified according to 
# http://techpubs.sgi.com/library/tpl/cgi-bin/getdoc.cgi/0640/bks/SGI_Developer/books/MProF90_PG/sgi_html/apb.html
# The two formulas (the one from the URL above and the one from lib.py) seem to
# agree only up to 4 deminal spaces, meaning to 1/3 of an arcsec.
dist = [173.011470, 142.995674, 119.905476, 115.707343, 156.220963,  84.231585,
        131.262275,  20.922682,  74.699713,  43.990162,  59.929432,  58.330367,  
        26.981633, 102.269873,  51.655353, 161.378036, 104.306836,  69.391945,  
        96.459234,  34.673781,  61.068864,  46.893570, 129.703251, 105.471314,  
        98.034748,  82.982402, 114.243931,  23.404158, 138.302827, 162.434956,  
        62.602563, 126.700180, 104.534651, 113.093854,  65.327320,  79.894500,  
        67.676184, 135.404569,  84.098437,  93.327769,  68.644661,  74.918530, 
        141.999974, 151.077623,  58.187654]
# These come from PS-MOPS and are based off of HI (UT offset -10)
mjds = [55015.061930981828, 55136.849253146684, 54924.207037036969, 
        55021.077138484106, 55187.069613991829, 55113.335028170186, 
        54893.086759953338, 55071.144613572716, 54913.534992714965, 
        54943.017330454306]
nn = [55014, 55135, 54923, 55020, 55186, 55112, 54892, 55070, 54912, 54942]
offset = -10.



class TestLib(unittest.TestCase):
    def setUp(self):
        return
    
    def testSphericalDistance(self):
        mydist = []
        n = len(ra)
        for i in range(n-1):
            for j in range(i+1, n, 1):
                mydist.append(lib.sphericalDistance((ra[i], dec[i]),
                                                    (ra[j], dec[j])))
        self.failUnlessEqual(len(dist), len(mydist), 'different number of distances.')
        
        # Test the total distances.
        for i in range(len(dist)):
            self.failUnlessAlmostEqual(dist[i], mydist[i][-1], 4, 'distances differ (%.04f /= %.04f).' %(dist[i], mydist[i][-1]))
        return
    
    def testMjdToNightNumber(self):
        for i in range(len(mjds)):
            mynn = lib.mjdToNightNumber(mjds[i], offset)
            self.failUnlessEqual(nn[i], mynn, 'night numbers differ: %d /= %d' %(nn[i], mynn))
        return
    



def sanityCheck():
    return



if(__name__ == '__main__'):
    sanityCheck()        
    unittest.main()





















'''
Created on Feb 18, 2015

@author: cmccully
'''
from pysao import ds9
from glob import glob
import pyfits
import numpy as np
from astropy.wcs import WCS
ra = 69.561622
dec = -12.287507


d = ds9.ds9(start=True)

fs = glob('*e90.fits')

fs = np.array(fs)

mjds = []

for f in fs:
    mjds.append(float(pyfits.getval(f, 'MJD-OBS')))

mjds = np.array(mjds)

fs = fs[np.argsort(mjds)]
for i in range(len(fs)):
    if i != 0:
        d.set('frame new')
    d.set('frame frameno %i'%(i+1))
    d.set('file %s'%fs[i])
    hdu = pyfits.open(fs[i])
    
    low = np.median(hdu[0].data) - 100
    w = WCS(fs[i])
    x,y = w.wcs_world2pix([ra], [dec], 0)
    high = 200 + hdu[0].data[int(y),int(x)]
    hdu.close()
    d.set('scale limits %f %f'%(low, high))

d.set('frame lock wcs')
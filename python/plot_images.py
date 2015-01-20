import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy import wcs
from wcsaxes import WCS
import aplpy

#ra = 04:38:14.9000 
ra = (4.+38./60.+14.9/60./60.)/24*360 
#dec = -12:17:14.40
dec = -12. - (17./60 +14.4/60./60.)
print 'ra',ra,'dec',dec
plt.rcParams['figure.figsize']=(10,8)
#fits_file_name = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0123-e90.fits'
#fits_file_name = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0124-e90.fits'
fits_file_name = '../data/HE0435_LCOGT/lsc1m009-fl03-20141201-0115-e90.fits'

hdulist=fits.open(fits_file_name)
hdulist.info()
bigw=WCS(hdulist[0].header)
image_data = hdulist[0].data

#w = WCS('image.fits')
lon, lat = bigw.all_pix2world(30, 40, 0)
print(lon, lat)
x,y = bigw.wcs_world2pix(ra,dec,1)
print x,y
lon, lat = bigw.all_pix2world(x, y, 0)
print(lon, lat)
#exit()

print(type(image_data))
print(image_data.shape)

#ax=plt.subplot(221)
NBINS = 1000
histogram = plt.hist(image_data.flat, log=True, bins=NBINS)
#plt.subplot(222)
plt.figure()
plt.imshow((image_data), cmap='gray', vmin=1400.,vmax=1800., interpolation='none')
#plt.imshow((image_data), cmap='gray', vmin=1400.,vmax=1800.)
plt.figure()
plt.imshow((image_data), cmap='gray', vmin=1.5*1400.,vmax=2.5*1800., interpolation='none')
plt.colorbar()
plt.xlim(2011-100,2011+100)
plt.ylim(2048+100,2048-100)

#plt.figure()
gc = aplpy.FITSFigure(fits_file_name)
gc.show_grayscale(vmin=np.percentile(image_data,5),vmax=np.percentile(image_data,98.5))
plt.show()

plt.colorbar()
#plt.subplot(223)
plt.imshow((image_data), cmap='gray', vmin=0.,vmax=8000.)
#plt.subplot(224)
plt.imshow((image_data), cmap='gray', vmin=0.,vmax=16000.)
#plt.xlim(1900.,2000.)
#plt.ylim(1900.,2000.)
plt.colorbar()
#plt.show()

plt.show()


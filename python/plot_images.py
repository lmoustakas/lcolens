import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy import wcs
from wcsaxes import WCS
import aplpy
from astropy.modeling import models, fitting
from astropy.modeling.models import Gaussian2D
import scipy.optimize as opt
from AnalysisLibrary import *


#ra = 04:38:14.9000 
ra = (4.+38./60.+14.9/60./60.)/24*360 
#dec = -12:17:14.40
dec = -12. - (17./60 +14.4/60./60.)
print 'ra',ra,'dec',dec
plt.rcParams['figure.figsize']=(10,8)
#fits_file_name = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0123-e90.fits'
#fits_file_name = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0124-e90.fits'
#fits_file_name = '../data/HE0435_LCOGT/lsc1m009-fl03-20141201-0115-e90.fits'
#fits_file_name = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0123-e90.fits'
#fits_file_name = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0123-e90_wcs_corrected.fits'
#fits_file_name = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0124-e90_wcs_corrected.fits'
fits_file_name = '../data/HE0435_LCOGT/lsc1m009-fl03-20141201-0115-e90_wcs_corrected.fits'

a1 = AnalysisLibrary(fits_file_name)

hdulist=fits.open(fits_file_name)
print 'hdulist.info()'
hdulist.info()
image_data = hdulist[0].data

print 'test bigw'
bigw=WCS(hdulist[0].header)
lon, lat = bigw.all_pix2world(30, 40, 0)
print(lon, lat)
x,y = bigw.wcs_world2pix(ra,dec,1)
print x,y
lon, lat = bigw.all_pix2world(x, y, 0)
print(lon, lat)

print 'astropy wcs'
astropy_wcs=wcs.WCS(hdulist[0].header)
#astropy_wcs=WCS(hdulist[0].header)
lon, lat = bigw.all_pix2world(30, 40, 0)
print(lon, lat)
x,y = bigw.wcs_world2pix(ra,dec,1)
print x,y
lon, lat = bigw.all_pix2world(x, y, 0)
print(lon, lat)
exit()

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

#plt.figure()
gc2 = aplpy.FITSFigure(fits_file_name)
gc2.recenter(ra, dec, radius=0.003) 
gc2.show_grayscale(vmin=np.percentile(image_data,5),vmax=np.percentile(image_data,98.5))
plt.figure()
x,y = bigw.wcs_world2pix(ra,dec,1)
plt.imshow(image_data, cmap='gray', vmin=np.percentile(image_data,5),vmax=np.percentile(image_data,98.5), interpolation='none')
plt.xlim(x-20,x+20)
plt.ylim(y-20,y+20)

plt.figure()
zoom_data = image_data[y-20:y+20,x-20:x+20]
#zoom_data = np.rot90(zoom_data)
zoom_data = np.flipud(zoom_data)
print zoom_data
print len(zoom_data)
#plt.imshow((image_data), cmap='gray', vmin=1.5*1400.,vmax=2.5*1800., interpolation='none')
plt.imshow((zoom_data), cmap='gray', interpolation='none')
plt.show()
print np.argmax(zoom_data)
#print np.argmax(zoom_data,axis=0)
#print np.argmax(zoom_data,axis=1)
print np.unravel_index(zoom_data.argmax(), zoom_data.shape)
x_max, y_max=np.unravel_index(zoom_data.argmax(), zoom_data.shape)
zoom_data2 = zoom_data[x_max-10:x_max+10,y_max-10:y_max+10]
#zoom_data2-=1500.
x_max2, y_max2=np.unravel_index(zoom_data2.argmax(), zoom_data2.shape)
print zoom_data2[x_max2,y_max2], zoom_data[y_max2,x_max2]
plt.figure()
plt.imshow((zoom_data2), cmap='gray', interpolation='none')
plt.colorbar()
#plt.show()
r=[]
a=[]
print 'flattening around peak'
print zoom_data2.shape
for ix in range(-9,9):
  for iy in range(-9,9):
    #print x_max2+ix,y_max2+iy
    rad=np.sqrt(ix**2+iy**2)
    val = zoom_data2[x_max2+ix,y_max2+iy]
    r.append(rad)
    a.append(val)
plt.figure()
plt.plot(r,a,'.')
plt.grid(True)
#plt.xlim(x-10,x+10)
#plt.ylim(y-10,y+10)

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
  xo = float(xo)
  yo = float(yo)    
  a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
  b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
  c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
  g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
			  + c*((y-yo)**2)))
  if(theta<=0.or theta>2.*np.pi): g+=1.e99
  #print np.array(g)
  #print g
  return g.ravel()
  #return np.array(g)
print 'fitting'
# Fit the data using astropy.modeling
#OffsetGaussian = Gaussian2D + Gaussian2D
#print OffsetGaussian
initial_guess = (zoom_data2[x_max2,y_max2],x_max2,y_max2,3.,3.,0.,1500.)
x=np.arange(0.,20)
y=np.arange(0.,20)
xg, yg = np.mgrid[:20, :20]
popt, pcov = opt.curve_fit(twoD_Gaussian, (xg, yg), zoom_data2.ravel(), p0=initial_guess)

print popt
fitted_data=twoD_Gaussian((xg,yg),*popt)
plt.figure()
plt.subplot(2,2,1)
plt.imshow(fitted_data.reshape(20,20), cmap='gray', interpolation='none')
plt.colorbar()
plt.subplot(2,2,2)
plt.imshow(zoom_data2-fitted_data.reshape(20,20), cmap='gray', interpolation='none')
plt.colorbar()
#plt.imshow(p_init(xg,yg), cmap='gray', interpolation='none')
plt.subplot(2,2,3)
x_rand=85
y_rand=135
plt.imshow(zoom_data[x_rand-10:x_rand+10,y_rand-10:y_rand+10], cmap='gray', interpolation='none')
plt.colorbar()

plt.show()


import numpy as np
from astropy.io import fits
from astropy import wcs

from wcsaxes import WCSAxes
import aplpy
import pylab as plt
import scipy.optimize as opt


class FITSmanager:
  def __init__(self,fits_file_name):
    self.fits_file_name = fits_file_name
    self.hdulist = fits.open(fits_file_name)
    #self.hdulist.info()
    self.image_data = self.hdulist[0].data
    self.bigw=wcs.WCS(self.hdulist[0].header)

  def histogram_image_values(self,NBINS=5000):
    self.hist, self.bin_edges = np.histogram(self.image_data.flat, bins=NBINS, range=(0.,200000.))
    self.BW = self.bin_edges[1] - self.bin_edges[0]


  def plot_image_values(self,NBINS=5000):
    #histogram = plt.hist(self.image_data.flat, log=True, bins=NBINS)
    self.histogram_image_values(NBINS=NBINS)
    plt.plot(self.bin_edges[:-1]+self.BW/2.,self.hist+1.e-1, label=self.fits_file_name)
    plt.xlabel('Image Value')
    plt.ylabel('Counts')
    
  def plot_image(self, ra_center=None, dec_center=None, rad=None, full_range=False):
    gc = aplpy.FITSFigure(self.fits_file_name)
    if(ra_center==None and dec_center==None and rad==None and full_range==True):
      print 'It is True, it is True, full_range is True'
      gc.show_grayscale()
    if(ra_center==None and dec_center==None and rad==None and full_range==False):
      gc.show_grayscale(vmin=np.percentile(self.image_data,5),vmax=np.percentile(self.image_data,98.5))
    if(ra_center!=None and dec_center!=None and rad!=None and full_range==False):
      gc.recenter(ra_center, dec_center, radius=rad) 
      gc.show_grayscale(vmin=np.percentile(self.image_data,5),vmax=np.percentile(self.image_data,98.5))
    gc.add_colorbar()
    xlabel='Right Ascension (J2000)'
    ylabel='Declination (J2000)'
    
  def image_piece(self,ra, dec, pixels):
    if(pixels%2!=0):
      pixels+=1
    x,y = self.bigw.wcs_world2pix(ra,dec,1)
    zoom_data = self.image_data[y-pixels/2:y+pixels/2,x-pixels/2:x+pixels/2]
    return zoom_data
    #plt.figure()

def twoD_Moffat((x, y), amplitude, alpha, beta, xo, yo, offset):
  xo = float(xo)
  yo = float(yo)    
  a = (beta-1.)/(np.pi*alpha**2)
  m = offset + amplitude*( 1. + ((x-xo)**2 + (y-yo)**2) / (2.*alpha**2))**(-beta)
  #print np.array(g)
  #print g
  print 'in Moffat 2D', offset, amplitude, a, amplitude*a
  if(alpha<0.): m+=1.e9
  #print offset
  return m.ravel()
  #return np.array(g)    

def twoD_Moffat_proj(x, amplitude, alpha, beta, xo, yo, offset):
  a = (beta-1.)/(np.pi*alpha**2)
  m = offset + amplitude * a * ( 1. + (x**2) / (2.*alpha**2))**(-beta)
  #print np.array(g)
  #print g
  return m.ravel()
  #return np.array(g)    
  
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

class SourceImage:
  def __init__(self,AnLib, ra,dec, pixels):
    self.image = AnLib.image_piece(ra, dec, pixels)
    self.ra = ra
    self.dec = dec
    self.pixels = pixels

  def fit_moffat(self, verbose=False):
    if(verbose): print 'fitting'
    self.x_max, self.y_max = np.unravel_index(self.image.argmax(), self.image.shape)
    background_count_guess = np.min(self.image)
    #initial_guess = (image[x_max,y_max],x_max,y_max,3.,3.,0.,1500.)
    #initial_guess_simple = (image[x_max,y_max],x_max,y_max,3.,1500.)
    guess_alpha = 3.
    guess_beta = 2.
    self.initial_guess_moffat = (self.image[self.x_max,self.y_max], guess_alpha, guess_beta, self.x_max, self.y_max, background_count_guess)
    if(verbose): print 'initial guess', self.initial_guess_moffat
    if(verbose): print len(self.image)
    x=np.arange(0.,len(self.image))
    y=np.arange(0.,len(self.image))
    self.xg, self.yg = np.mgrid[:len(self.image), :len(self.image)]
    #popt, pcov = opt.curve_fit(twoD_Gaussian, (xg, yg), image.ravel(), p0=initial_guess)
    #popt, pcov = opt.curve_fit(twoD_Gaussian_simple, (xg, yg), image.ravel(), p0=initial_guess_simple)
    popt, pcov = opt.curve_fit(twoD_Moffat, (self.xg, self.yg), self.image.ravel(), p0=self.initial_guess_moffat)
    if(verbose): print popt
    
    self.moffat_fit_image = twoD_Moffat((self.xg, self.yg), *popt).reshape(len(self.image),len(self.image))

    self.rad, self.amp = radial_profile(self.image, popt[3], popt[4])
    self.moffat_parms = popt
    self.moffat_fwhm = popt[1]*2.*np.sqrt(2.**(1/popt[2])-1.)
    self.moffat_chi = (self.amp-twoD_Moffat_proj(self.rad, *self.moffat_parms))/np.sqrt(self.amp)
    #self.moffat_chi = (self.amp-twoD_Moffat_proj(self.rad, *self.moffat_parms))/np.sqrt(twoD_Moffat_proj(self.rad, *self.moffat_parms))
    return popt

  def plot_moffat_residual_2D(self):
    plt.imshow(self.image-self.moffat_fit_image, cmap='gray', interpolation='none')
    plt.title('ra: %1.4f dec: %1.4f'%(self.ra,self.dec))
    plt.colorbar()

  def plot_radial_profiles(self):
    plt.plot(self.rad,self.amp,'b.')
    txt = 'amplitude: %1.2e\n'%self.moffat_parms[0]
    txt+= 'alpha:        %1.2f\n'%self.moffat_parms[1]
    txt+= 'beta:          %1.2f\n'%self.moffat_parms[2]
    txt+= 'FWHM:          %1.2f\n'%self.moffat_fwhm
    txt+= 'offset:       %1.2e'%self.moffat_parms[5]
    #plt.plot(self.rad,twoD_Moffat_proj(self.rad, *self.moffat_parms),'r-',label=txt)
    x = np.arange(0.,self.pixels,0.01)
    plt.plot(x,twoD_Moffat_proj(x, *self.moffat_parms),'r-',label=txt)
    plt.legend(loc=1)
    plt.title('ra: %1.4f dec: %1.4f'%(self.ra,self.dec))
    plt.xlabel('Pixel Distance from Fitted Peak')
    plt.ylabel('Pixel Value')

  def plot_radial_profile_residuals(self):
    plt.plot(self.rad,self.moffat_chi,'b.')
    plt.title('ra: %1.4f dec: %1.4f'%(self.ra,self.dec))
    plt.ylim(-15.,15.)
    plt.xlabel('Pixel Distance from Fitted Peak')
    plt.ylabel(r'$\chi_p$', fontsize=20)

  def plot_moffat_chi(self):
    chisq = np.sum((self.moffat_chi)**2)/(len(self.rad)-6.)
    a,b,c = plt.hist(self.moffat_chi,bins=51, range=(-10,10), label=r'$\chi^2$=%1.2f'%chisq)
    #plt.legend(loc=1)
    label='$\chi^2$/d.o.f.=%1.2f'%chisq
    plt.text(2.,0.9*max(a), label, size=16)

    plt.title('ra: %1.4f dec: %1.4f'%(self.ra,self.dec))
    plt.xlim(-10.,10.)
    #plt.ylim(-15.,15.)
    #plt.xlabel('Pixel Distance from Fitted Peak')
    plt.xlabel(r'$\chi_p$')

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

def twoD_Gaussian_simple((x, y), amplitude, xo, yo, sigma, offset):
  xo = float(xo)
  yo = float(yo)    
  a = 1./(2*sigma**2)
  g = offset + amplitude*np.exp( - ((x-xo)**2 + (y-yo)**2) / (2.*sigma**2))
  #print np.array(g)
  #print g
  return g.ravel()
  #return np.array(g)  
  


def fit_guassian(image, verbose=False):

  if(verbose): print 'fitting'
  x_max, y_max=np.unravel_index(image.argmax(), image.shape)

  #initial_guess = (image[x_max,y_max],x_max,y_max,3.,3.,0.,1500.)
  #initial_guess_simple = (image[x_max,y_max],x_max,y_max,3.,1500.)
  initial_guess_moffat = (image[x_max,y_max],3.,2.,x_max,y_max,1500.)
  if(verbose):print len(image)
  x=np.arange(0.,len(image))
  y=np.arange(0.,len(image))
  xg, yg = np.mgrid[:len(image), :len(image)]
  #popt, pcov = opt.curve_fit(twoD_Gaussian, (xg, yg), image.ravel(), p0=initial_guess)
  #popt, pcov = opt.curve_fit(twoD_Gaussian_simple, (xg, yg), image.ravel(), p0=initial_guess_simple)
  popt, pcov = opt.curve_fit(twoD_Moffat, (xg, yg), image.ravel(), p0=initial_guess_moffat)
  if(verbose):print popt
  return popt

def radial_profile(image, x0,y0):
  r=[]
  amp=[]
  for i in range(0,len(image)):
      for j in range(0,len(image)):
	r.append(np.sqrt((float(i)-x0)**2 + (float(j)-y0)**2))
	amp.append(image[i][j])
  return np.array(r), np.array(amp)


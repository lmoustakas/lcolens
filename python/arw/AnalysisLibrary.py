import numpy as np
import datetime
import os
from astropy.io import fits
from astropy.io import ascii
from astropy import wcs

from wcsaxes import WCSAxes
import aplpy
import pylab as plt
import scipy.optimize as opt
from scipy.optimize import minimize


def time_order_fits_files(data_dir):
	# GET ALL DATA FILE NAMES IN DATA DIRECTORY AND SORT THEM BY MJD
	print 'READING THE FITS FILE OBSERVATION TIMES FOR ORDERING PURPOSES'
	fnms = os.listdir(data_dir)
	print 'Number of files in directory:', len(fnms)
	fnms_filtered = []
	for f in fnms:
	  if(len(f.split('-'))>3):
	     fnms_filtered.append(f)
	fnms = fnms_filtered
	print 'Number of files after filtering for valid data file names:', len(fnms)

	print 'READING THE FITS FILE OBSERVATION TIMES FOR ORDERING PURPOSES'
	mjd_obs_list = []
	fnms_list=[]
	for i in range(len(fnms)):
	    if(i%50==0): print '%d of %d files read'%(i,len(fnms))
	    d=fnms[i].split('-')
	    #print i,fnms[i]
	    fnms_list.append(fnms[i])
	    fits_file_name = data_dir + fnms[i]
	    hdulist = fits.open(fits_file_name)
	    mjd_obs_list.append(float(hdulist[0].header['MJD-OBS']))
	    hdulist.close()

	print 'SORTING THE FILE NAME LIST IN TIME ORDER'

	fnms    = [x for (y,x) in sorted(zip(mjd_obs_list, fnms_list))]
	fnm_mjd = [y for (y,x) in sorted(zip(mjd_obs_list, fnms_list))]
	print 'Number of files after sorting:', len(fnms)
	ascii.write([fnms, fnm_mjd], 'time_ordered_file_list.dat', names=['filename', 'mjd'] )

	return fnms

################################################################################################################
################################################################################################################
################################################################################################################

class FITSmanager:
  def __init__(self,fits_file_name):
    self.fits_file_name = fits_file_name
    self.hdulist = fits.open(fits_file_name)
    #self.hdulist.info()
    self.image_data = self.hdulist[0].data
    self.bigw=wcs.WCS(self.hdulist[0].header)

  def histogram_image_values(self,NBINS=5000, rng_max=200000.):
    self.hist, self.bin_edges = np.histogram(self.image_data.flat, bins=NBINS, range=(0.,rng_max))
    self.BW = self.bin_edges[1] - self.bin_edges[0]


  def plot_image_values(self,NBINS=5000, rng_max=200000.):
    #histogram = plt.hist(self.image_data.flat, log=True, bins=NBINS)
    self.histogram_image_values(NBINS=NBINS, rng_max=rng_max)
    plt.plot(self.bin_edges[:-1]+self.BW/2.,self.hist+1.e-1, label=self.fits_file_name)
    plt.xlabel('Image Value')
    plt.ylabel('Counts')
    
  def plot_image(self, ra_center=None, dec_center=None, rad=None, full_range=False):
    gc = aplpy.FITSFigure(self.fits_file_name)
    if(ra_center==None and dec_center==None and rad==None and full_range==True):
      print 'It is True, it is True, full_range is True'
      gc.show_colorscale(vmin=np.percentile(self.image_data,0.1),vmax=np.percentile(self.image_data,99.9))
      plt.title('Full Range')
      gc.add_colorbar()
      xlabel='Right Ascension (J2000)'
      ylabel='Declination (J2000)'
      return
    if(ra_center==None and dec_center==None and rad==None and full_range==False):
      gc.show_grayscale(vmin=np.percentile(self.image_data,5),vmax=np.percentile(self.image_data,98.5))
    if(ra_center!=None and dec_center!=None and rad!=None and full_range==False):
      gc.recenter(ra_center, dec_center, radius=rad) 
      gc.show_grayscale(vmin=np.percentile(self.image_data,5),vmax=np.percentile(self.image_data,98.5))
    gc.add_colorbar()
    xlabel='Right Ascension (J2000)'
    ylabel='Declination (J2000)'
    
  def image_piece(self,ra, dec, pixels):
    x,y = self.bigw.wcs_world2pix(ra,dec,1)
    if(pixels%2==0):
	    #print 'even'
	    zoom_data = self.image_data[y-pixels/2:y+pixels/2,x-pixels/2:x+pixels/2]
    if(pixels%2==1):
	    #print 'odd'
	    zoom_data = self.image_data[y-(pixels-1)/2:y+(pixels-1)/2+1,x-(pixels-1)/2:x+(pixels-1)/2+1]
    return zoom_data
    #plt.figure()

  def get_exposure_time(self):
	  t1 = self.hdulist[0].header['UTSTART']
	  t2 = self.hdulist[0].header['UTSTOP']
	  #print t1+'000',t2
	  self.start_time = datetime.datetime.strptime(t1,"%H:%M:%S.%f")
	  self.end_time = datetime.datetime.strptime(t2,"%H:%M:%S.%f")
	  self.exposure_time = (self.end_time - self.start_time).seconds + (self.end_time - self.start_time).microseconds/1.e6
	  return self.exposure_time


def twoD_Moffat((x, y), amplitude, alpha, beta, xo, yo, offset):
  xo = float(xo)
  yo = float(yo)    
  a = (beta-1.)/(np.pi*alpha**2)
  m = offset + amplitude*( 1. + ((x-xo)**2 + (y-yo)**2) / (2.*alpha**2))**(-beta)
  #print np.array(g)
  #print g
  #print 'in Moffat 2D', offset, amplitude, a, amplitude*a
  if(alpha<0.): m+=1.e9
  #print offset
  return m.ravel()
  #return np.array(g)    

def twoD_Moffat_proj(x, amplitude, alpha, beta, xo, yo, offset):
  a = (beta-1.)/(np.pi*alpha**2)
  m = offset + amplitude * ( 1. + (x**2) / (2.*alpha**2))**(-beta)
  #print np.array(g)
  #print g
  return m.ravel()
  #return np.array(g)    
  
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

class SourceImage:
  def __init__(self,FM, ra,dec, pixels):
    self.FM = FM
    self.image = FM.image_piece(ra, dec, pixels)
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




################################################################
################################################################
################################################################

def estimate_background_and_read_noise(fits_file_name, display=False):
	# READ FITS FILE
	FM = FITSmanager(fits_file_name)
	# PRODUCE A HISTOGRAM OF THE IMAGE COUNTS IN THE LOW END
	FM.histogram_image_values(NBINS=10000, rng_max=10000.)
	# ESTIMATE THE BACKGROUND COUNTS AS THE MOST COMMON VALUE OF THE IMAGE
	N_bkg = FM.bin_edges[np.argmax(FM.hist)]
	# GET THE LOWER POINT AT WHICH THE COUNTS ARE REDUCED BY A FACTOR OF EXP(-0.5) 
	sig_N_bkg_index = np.argmin((FM.hist[:np.argmax(FM.hist)]-np.exp(-0.5)*np.max(FM.hist))**2) 
	# ESTIMATE THE SIGMA OF THE DISTRIBUTION BASED ON THE LOCATION OF THIS POINT 
	sig_N_bkg = FM.bin_edges[np.argmax(FM.hist)] -  FM.bin_edges[sig_N_bkg_index]
	# THE READ NOISE IS THE DISCREPANCY (in quadrature) BETWEEN sqrt(N_bkg) AND THE SIGMA OF THE DISTRIBUTION ESTIMATED ABOVE
	sigma_read = 0.
	if(sig_N_bkg**2-FM.bin_edges[np.argmax(FM.hist)] > 0.):
		sigma_read = np.sqrt(sig_N_bkg**2-FM.bin_edges[np.argmax(FM.hist)])
	print 'Pixel Value Background:\t%1.0f'%FM.bin_edges[np.argmax(FM.hist)]
	print 'Bkg Standard Deviation:\t%1.2f'%sig_N_bkg
	print 'sqrt(N_bkg)\t\t%1.2f'%np.sqrt(FM.bin_edges[np.argmax(FM.hist)])
	print 'Read noise\t\t%1.2f'%sigma_read
	if(display):

		#subplot(211)
		plt.plot(FM.bin_edges[:-1],FM.hist, lw=2)
		plt.plot(FM.bin_edges[:np.argmax(FM.hist)],FM.hist[:np.argmax(FM.hist)], '--', lw=2, color='orange')
		plt.plot([FM.bin_edges[np.argmax(FM.hist)],FM.bin_edges[np.argmax(FM.hist)]],[0.,max(FM.hist)],'r--', lw=2)
		plt.plot([FM.bin_edges[sig_N_bkg_index],FM.bin_edges[sig_N_bkg_index]],[0.,np.exp(-0.5)*max(FM.hist)],'g--', lw=2)
		plt.xlim(0., 3.*FM.bin_edges[np.argmax(FM.hist)])
		plt.ylabel('Counts')
		plt.xlabel('Pixel Value')
		#print k, FM.bin_edges[np.argmax(FM.hist)]
		plt.show()
	return N_bkg, sigma_read

def make_background_counts_and_read_noise_table():
  dirc = '/data2/romerowo/lcogt_data/he045-1223_wcs_corrected/'
  filename_table = ascii.read('time_ordered_file_list.dat')
  N_bkg = []
  sigma_read = []
  count=0
  bad_observation_list=[364,365, 386]
  #bad_observation_list=[-1]
  print len(filename_table['filename'])
  for k in range(0,len(filename_table['filename'])):
    FM = FITSmanager(dirc+filename_table['filename'][k])
    #FM.plot_image_values(NBINS=10000, rng_max=10000)
    print '\n'
    print '%d of %d: %s'%(k, len(filename_table['filename']), filename_table['filename'][k])
    if(k in bad_observation_list):
	print ' The file %s has been tagged as bad data'%filename_table['filename'][k]
	N_bkg.append(-1)
	sigma_read.append(-1)
    if(k not in bad_observation_list):
	N_bkg_val, sigma_read_val = estimate_background_and_read_noise(dirc+filename_table['filename'][k])
	N_bkg.append(N_bkg_val)
	sigma_read.append(sigma_read_val)
	count += 1

  ascii.write([filename_table['filename'], filename_table['mjd'], N_bkg, sigma_read], 'bkg_counts_and_read_noise.tbl', names=['filename', 'mjd', 'N_bkg', 'read_noise'])


def gaussian((x, y), amp, sigx, sigy, ang, x0, y0):
    # Precompute the sin and cos
    spa = np.sin(ang / 180.0 * np.pi)
    cpa = np.cos(ang / 180.0 * np.pi)
    #print spa, cpa
    # Define xprime coordinates in the rotated frame for convenience
    # This uses a standard rotation matrix
    xp = (x - x0) * cpa - (y - y0) * spa
    yp = (x - x0) * spa + (y - y0) * cpa
    # Defined r^2 (no need to take a square root)
    r2 = xp * xp / sigx / sigx + yp * yp / sigy / sigy
    return amp*np.exp(-r2*r2/2.)

def triple_gaussian((x, y), amp1, sigx1, sigy1, ang1, xo1, yo1, amp2, sigx2, sigy2, ang2, xo2, yo2, amp3, sigx3, sigy3, ang3, xo3, yo3):
	g1 = gaussian((x, y), amp1, sigx1, sigy1, ang1, xo1, yo1)
	g2 = gaussian((x, y), amp2, sigx2, sigy2, ang2, xo2, yo2)
	g3 = gaussian((x, y), amp3, sigx3, sigy3, ang3, xo3, yo3)
	return g1 + g2 + g3

def triple_gauss_chi2(parms, image, sig):
	xg, yg = np.mgrid[:len(image), :len(image)]
	model = triple_gaussian( (xg,yg,), *parms)
	chi = (image-model)/sig
	print np.sum(chi*chi), parms
	return np.sum(chi*chi)

def fitter(image, p0):
    sig = np.sqrt(image)
    results = minimize(triple_gauss_chi2, p0, args=(image, sig), method='Nelder-Mead')
    print 'finished running minimize'
    # return the best fit parameters
    print results
    return results.x
#get_background_counts()
#exit()


import numpy as np
import datetime
import os
from astropy.io import fits
from astropy.io import ascii
from astropy import wcs

#from wcsaxes import WCSAxes
import aplpy
import pylab as plt
import scipy.optimize as opt
from scipy.optimize import minimize
from scipy.optimize import leastsq
import matplotlib
from scipy import signal
import emcee
import triangle


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
    # MULTIPLY BY THE GAIN
    self.image_data *= float(self.hdulist[0].header['GAIN'])
    self.bigw=wcs.WCS(self.hdulist[0].header)

  def estimate_read_noise(self, display=0, out = 'out'):
	readnoise=[]
	for k in range(1000):
		#print self.image_data.shape
		x = np.random.randint(100, self.image_data.shape[0]-100)
		y = np.random.randint(100, self.image_data.shape[1]-100)
		r, d = self.bigw.wcs_pix2world(x,y, 1)
		#print 'x,y',x,y
		#print 'ra,dec',r,d
		img = self.image_piece(r, d, 31)
		#img/=float(self.hdulist[0].header['GAIN']) # remove gain correction. read noise is on electrons, not photons ?		
		count = 0
		
		while(np.max(img)>np.median(img)+5.*np.sqrt(np.median(img))):
			count+=1
	  		#print 'trial', count
			x = np.random.randint(100, self.image_data.shape[0]-100)
			y = np.random.randint(100, self.image_data.shape[1]-100)
			r, d = self.bigw.wcs_pix2world(x,y, 1)
			#print 'x,y',x,y
			#print 'ra,dec',r,d
			img = self.image_piece(r, d, 31)
		
		h,b = np.histogram(img.flat, bins=30)
		x=b[:-1]+(b[1]-b[0])/2.
		try:
			popt, pcov = opt.curve_fit(gauss_1d, (x), h, p0=[np.max(h), np.sqrt(np.median(x)), np.median(x)])
		except:
			continue
		if(popt[0]<=0. or popt[2]<0):
			continue
		if(popt[1]**2 - popt[2] <= 0.):
			continue
		rn = np.sqrt(popt[1]**2 - popt[2])
		if(rn == rn):
			readnoise.append(rn)
	 		#print k,'readnoise', popt[1]-np.sqrt(popt[0])
	 		#print 'median', np.median(img)
			#print popt
	self.readnoise = np.median(readnoise)
	print np.median(readnoise)
	if(display>=2):
		plt.figure(figsize=(9,7))
		plt.subplot(221)
		plt.imshow(img, interpolation='none')
		plt.colorbar()
		plt.title('Random Image Portion')
		plt.xlabel('pixel')
		plt.ylabel('pixel')
		plt.subplot(222)
		plt.step(x,h, lw=2)
		plt.plot(x,gauss_1d(x,*popt), 'r--', lw=2) 
		plt.title('Distribution of Counts\nw/ Gaussian Fit')
		plt.xlabel('CCD Counts (gain corrected)')
		plt.subplot(212)
		print 'median read noise', np.median(readnoise)
		a,b,c = plt.hist(readnoise, bins=int(np.max(readnoise)-np.min(readnoise)))
		plt.plot([np.median(readnoise), np.median(readnoise)],[0., 1.2*max(a)],'r-', lw=2, label='median')
		plt.legend(loc=1)
		plt.ylim(0.,1.2*max(a))
		plt.xlabel('Read Noise, CCD counts (gain corrected)')
		plt.title('Multiple Read Noise Estimates')
		plt.suptitle(self.fits_file_name.split('/')[-1], fontsize=20)
		plt.subplots_adjust(top=0.85, wspace=0.3, hspace=0.4)
		plt.savefig(out+'.png', dpi=50)
	#return self.readnoise
	

  def histogram_image_values(self,NBINS=5000, rng_max=200000.):
    self.hist, self.bin_edges = np.histogram(self.image_data.flat, bins=NBINS, range=(0.,rng_max))
    self.BW = self.bin_edges[1] - self.bin_edges[0]


  def plot_image_values(self,NBINS=5000, rng_max=200000.):
    #histogram = plt.hist(self.image_data.flat, log=True, bins=NBINS)
    self.histogram_image_values(NBINS=NBINS, rng_max=rng_max)
    plt.plot(self.bin_edges[:-1]+self.BW/2.,self.hist+1.e-1, label=self.fits_file_name)
    plt.xlabel('Image Value')
    plt.ylabel('Counts')

  def plot_image(self, ra, dec, Npx=31, out='out'):
    # GET THE QUASAR IMAGE
    x,y = self.bigw.wcs_world2pix(ra,dec,1)    
    plt.figure(figsize=(9.5,8))
    ax=plt.subplot(111)
    yg, xg = np.mgrid[y-Npx/2:y+Npx/2+1,x-Npx/2:x+Npx/2+1]
    d, r = self.bigw.wcs_pix2world(xg,yg, 1)
    #plt.pcolor(d,r,self.image_data[y-Npx/2:y+Npx/2+1,x-Npx/2:x+Npx/2+1], interpolation='none')
    plt.pcolor(d,r,self.image_data[y-Npx/2:y+Npx/2+1,x-Npx/2:x+Npx/2+1], cmap='jet')
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)
    plt.colorbar()
    plt.xlabel('Right Ascension, deg', fontsize=14)
    plt.ylabel('Declination, deg', fontsize=14)
    print self.fits_file_name
    plt.title(self.fits_file_name.split('/')[-1])
    plt.subplots_adjust(left=0.15, right=1.0)
    plt.axis('equal')
    plt.ylim(np.min(r), np.max(r))
    plt.xlim(np.min(d), np.max(d))
    plt.savefig(out+'.png', dpi=50)
    #plt.show()

    
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

def twoD_elliptical_Moffat((x, y), amplitude, alpha, beta, xo, yo, el, theta, offset):
  xo = float(xo)
  yo = float(yo)    
  sigma_x = alpha
  sigma_y = el*alpha
  a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
  b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
  c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)

  m = offset + amplitude*( 1. + ( a*(x-xo)**2 + 2*b*(x-xo)*(y-yo)+ c*(y-yo)**2 ) )**(-beta)
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
    #while( np.min(self.image)<0. ):
	#xm, ym = np.unravel_index(self.image.argmin(), self.image.shape)
	#self.image[xm,ym]=np.median(self.image)
    if(background_count_guess<0.): 
	print np.median(self.image)
	plt.figure()
	plt.imshow(self.image, interpolation='none')
	plt.colorbar()
	plt.savefig('wtf.png')
	exit()
    #initial_guess = (image[x_max,y_max],x_max,y_max,3.,3.,0.,1500.)
    #initial_guess_simple = (image[x_max,y_max],x_max,y_max,3.,1500.)
    guess_alpha = 3.
    guess_beta = 2.
    self.initial_guess_moffat = (self.image[self.x_max,self.y_max], guess_alpha, guess_beta, self.x_max, self.y_max, background_count_guess)
    print self.initial_guess_moffat
    if(verbose): print 'initial guess', self.initial_guess_moffat
    if(verbose): print len(self.image)
    x=np.arange(0.,len(self.image))
    y=np.arange(0.,len(self.image))
    self.xg, self.yg = np.mgrid[:len(self.image), :len(self.image)]
    #popt, pcov = opt.curve_fit(twoD_Gaussian, (xg, yg), image.ravel(), p0=initial_guess)
    #popt, pcov = opt.curve_fit(twoD_Gaussian_simple, (xg, yg), image.ravel(), p0=initial_guess_simple)
    self.fit_ok = True
    #try:	
    popt, pcov = opt.curve_fit(twoD_Moffat, (self.xg, self.yg), self.image.ravel(), p0=self.initial_guess_moffat)
    '''
    except:
      self.fit_ok = False
      #maskval = np.sqrt((3.*np.sqrt(N_bkg))**2 + sigma_read**2)
      #Hmasked  = np.ma.masked_where((obj.image-N_bkg)<maskval,obj.image-N_bkg)
      Hmasked2 = np.ma.masked_where( (self.xg-15)**2 + (self.yg-15)**2 > 200,self.image)
      self.image = Hmasked2.copy()
      #estimate = np.sum(Hmasked)
      #popt, pcov = opt.curve_fit(twoD_Moffat, (self.xg, self.yg), self.image.ravel(), p0=self.initial_guess_moffat)
      print 'Moffat Fit Failed'
      #plt.figure()
      #plt.subplot(221)
      #plt.imshow(self.image, interpolation='none')
      #plt.subplot(222)
      #plt.imshow(Hmasked2, interpolation='none')
      #plt.show()
      popt = self.initial_guess_moffat
    '''
    if(verbose): print popt
    
    self.moffat_fit_image = twoD_Moffat((self.xg, self.yg), *popt).reshape(len(self.image),len(self.image))

    self.rad, self.amp = radial_profile(self.image, popt[3], popt[4])
    self.moffat_parms = popt
    self.moffat_fwhm = popt[1]*2.*np.sqrt(2.**(1/popt[2])-1.)
    self.moffat_chi = (self.amp-twoD_Moffat_proj(self.rad, *self.moffat_parms))/np.sqrt(self.amp)
    #self.moffat_chi = (self.amp-twoD_Moffat_proj(self.rad, *self.moffat_parms))/np.sqrt(twoD_Moffat_proj(self.rad, *self.moffat_parms))
    return popt


  def fit_elliptical_moffat(self, verbose=False):
    if(verbose): print 'fitting'
    self.x_max, self.y_max = np.unravel_index(self.image.argmax(), self.image.shape)
    background_count_guess = np.min(self.image)
    #while( np.min(self.image)<0. ):
	#xm, ym = np.unravel_index(self.image.argmin(), self.image.shape)
	#self.image[xm,ym]=np.median(self.image)
    if(background_count_guess<0.): 
	print np.median(self.image)
	plt.figure()
	plt.imshow(self.image, interpolation='none')
	plt.colorbar()
	plt.savefig('wtf.png')
	exit()
    #initial_guess = (image[x_max,y_max],x_max,y_max,3.,3.,0.,1500.)
    #initial_guess_simple = (image[x_max,y_max],x_max,y_max,3.,1500.)
    guess_alpha = 3.
    guess_beta = 2.
    el=0.
    self.initial_guess_elliptical_moffat = (self.image[self.x_max,self.y_max], guess_alpha, guess_beta, self.x_max, self.y_max, 1., 0., background_count_guess)
    print self.initial_guess_elliptical_moffat
    if(verbose): print 'initial guess', self.initial_guess_elliptical_moffat
    if(verbose): print len(self.image)
    x=np.arange(0.,len(self.image))
    y=np.arange(0.,len(self.image))
    self.xg, self.yg = np.mgrid[:len(self.image), :len(self.image)]
    #popt, pcov = opt.curve_fit(twoD_Gaussian, (xg, yg), image.ravel(), p0=initial_guess)
    #popt, pcov = opt.curve_fit(twoD_Gaussian_simple, (xg, yg), image.ravel(), p0=initial_guess_simple)
    self.fit_ok = True
    #try:	
    popt, pcov = opt.curve_fit(twoD_elliptical_Moffat, (self.xg, self.yg), self.image.ravel(), p0=self.initial_guess_elliptical_moffat)
    '''
    except:
      self.fit_ok = False
      #maskval = np.sqrt((3.*np.sqrt(N_bkg))**2 + sigma_read**2)
      #Hmasked  = np.ma.masked_where((obj.image-N_bkg)<maskval,obj.image-N_bkg)
      Hmasked2 = np.ma.masked_where( (self.xg-15)**2 + (self.yg-15)**2 > 200,self.image)
      self.image = Hmasked2.copy()
      #estimate = np.sum(Hmasked)
      #popt, pcov = opt.curve_fit(twoD_Moffat, (self.xg, self.yg), self.image.ravel(), p0=self.initial_guess_moffat)
      print 'Moffat Fit Failed'
      #plt.figure()
      #plt.subplot(221)
      #plt.imshow(self.image, interpolation='none')
      #plt.subplot(222)
      #plt.imshow(Hmasked2, interpolation='none')
      #plt.show()
      popt = self.initial_guess_moffat
    '''
    if(verbose): print popt
    
    self.elliptical_moffat_fit_image = twoD_elliptical_Moffat((self.xg, self.yg), *popt).reshape(len(self.image),len(self.image))
    self.rad, self.amp = radial_profile(self.image, popt[3], popt[4])
    self.elliptical_moffat_parms = popt
    self.elliptical_moffat_fwhm = popt[1]*2.*np.sqrt(2.**(1/popt[2])-1.)
    self.elliptical_moffat_chi = self.amp
    #self.elliptical_moffat_chi = (self.amp-twoD_Moffat_proj(self.rad, *self.elliptical_moffat_parms))/np.sqrt(self.amp)
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
  def quad_image_model(self, x0, y0, amp0, amp1, amp2, amp3, alpha, beta, N_bkg, N_pix, flip):
    #print len(theta), N_bkg, N_pix

    # NEED A GOOD MODEL HERE
    scale = 0.99*0.467/self.FM.hdulist[0].header['PIXSCALE']
    '''
    x1 = 14.3*scale + x0
    y1 = 13.2*scale + y0
    x2 = 15.6*scale + x0
    y2 = 18.4*scale + y0
    x3 = 13.1*scale + x0
    y3 = 16.4*scale + y0
    x4 = 17.6*scale + x0
    y4 = 15.3*scale + y0
    '''
    x1 = (14.3 + x0)*scale
    y1 = (13.2 + y0)*scale
    x2 = (15.6 + x0)*scale
    y2 = (18.4 + y0)*scale
    x3 = (13.1 + x0)*scale
    y3 = (16.4 + y0)*scale
    x4 = (17.6 + x0)*scale
    y4 = (15.3 + y0)*scale
    

    xg, yg = np.mgrid[:N_pix,:N_pix]

    #print x0,y0, amp0
    p0  =  twoD_Moffat((xg, yg), amp0, alpha, beta, x1, y1, 0)
    p1  =  twoD_Moffat((xg, yg), amp1, alpha, beta, x2, y2, 0)
    p2  =  twoD_Moffat((xg, yg), amp2, alpha, beta, x3, y3, 0)
    p3  =  twoD_Moffat((xg, yg), amp3, alpha, beta, x4, y4, 0)
    model = (p0+p1+p2+p3).reshape(N_pix, N_pix)
    if(flip):
    	model = np.fliplr(model)
    model+=N_bkg
    return model

  def moffat_chi_vals(self,theta,N_pix, flip):
    x0,y0,amp0,amp1,amp2,amp3, alpha, beta, N_bkg = theta
    model = self.quad_image_model(x0,y0,amp0,amp1,amp2,amp3, alpha, beta, N_bkg, N_pix, flip)
    chi = (self.image - model)/np.sqrt(self.image+self.FM.readnoise**2)
    #print np.sum(chi*chi)
    #print x0,y0, N_pix/2
    # NO NEGATIVE AMPLITUDES ALLOWED!
    if(amp0<0 or amp1<0 or amp2<0 or amp3<0):
	return np.inf
    # THE IMAGE CORRECTION HAS TO BE WITHIN THE BOUNDS OF THE IMAGE!
    if(x0>N_pix/2 or x0<-N_pix/2 or y0>N_pix/2 or y0<-N_pix/2):
	return np.inf
    return chi

  def moffat_chi_sq(self, theta, N_pix, flip):
	chisq = np.sum((self.moffat_chi_vals(theta, N_pix, flip))**2)
	#print '%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e'%(chisq, theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], theta[6], theta[7], theta[8]) 
	return chisq

  def moffat_chi_vals2(self,theta,alpha, beta,N_pix, flip):
    x0,y0,amp0,amp1,amp2,amp3, N_bkg = theta
    model = self.quad_image_model(x0,y0,amp0,amp1,amp2,amp3, alpha, beta, N_bkg, N_pix, flip)

    chi = (self.image - model)/np.sqrt(self.image+self.FM.readnoise**2)
    #print 'self.FM.readnoise', self.FM.readnoise
    #print np.sum(chi*chi)
    return chi

  def moffat_chi_sq2(self, theta, alpha, beta, N_pix, flip):
	return np.sum((self.moffat_chi_vals2(theta, alpha, beta, N_pix, flip))**2)


###############################################################################
## END CLASS DEFINITIONS ######################################################
###############################################################################


def moffat_analytical_integral(amplitude, alpha, beta):
  return amplitude * 2. * np.pi * alpha**2 / (beta-1.)


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
		#plt.show()
	return N_bkg, sigma_read

def moffat_fwhm(alpha,beta):
	return alpha*2.*np.sqrt(2.**(1/beta)-1.)
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


def gauss_1d(x, amp, sig, x0):
	return amp*np.exp(-(x-x0)**2/2./sig**2)

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

####################################################################
## PHOTOMETRY FUNCTIONS ############################################
####################################################################

def estimate_total_light(obj, N_bkg, sigma_read, display=0, out='out'):
  obj.fit_moffat()
  #print obj.moffat_parms
  alpha = obj.moffat_parms[1]
  beta = obj.moffat_parms[2]
  #print 'a,b',alpha,beta
  # DO SOME COOKIE CUT OUT STUFF LATER, FOR NOW ESTIMATE BY THE NUMBER OF COUNTS DUE TO THE SIGNAL
  estimate = np.sum(obj.moffat_fit_image - obj.moffat_parms[5]) # magnitude based on value of fitted counts
  uncertainty = np.sqrt(np.sum(obj.moffat_fit_image)+sigma_read**2)
  # percent uncertainty is on the level of light.
  #print '\t',N_bkg, obj.moffat_parms[5]
  #maskval = np.sqrt((3.*np.sqrt(N_bkg))**2 + sigma_read**2)
  #maskval = np.sqrt((5.*np.sqrt(obj.moffat_parms[5]))**2 + sigma_read**2)
  #Hmasked  = np.ma.masked_where((obj.image-obj.moffat_parms[5])<maskval,obj.image-obj.moffat_parms[5])
  #Hmasked2 = np.ma.masked_where((obj.image-obj.moffat_parms[5])<maskval,obj.image)
  #estimate = np.sum(Hmasked)
  #frac_uncertainty = np.sqrt(np.sum(Hmasked2)+sigma_read**2) / np.sum(Hmasked)
  print '\timage peak',np.max(obj.image)
  print '\timage bkg',np.min(obj.image)
  #print '\tmasked image sum counts',np.sum(Hmasked2)
  #print '\tmasked bkg subtractected image sum counts',np.sum(Hmasked)
  #print '\tfrac_uncertainty',frac_uncertainty
  print '\tmoffat_integrated count',estimate
  print '\tmoffat estimated uncertainty',uncertainty
  print '\tmoffat estimated frac uncertainty',uncertainty/estimate
  
  #maskval = np.sqrt((3.*np.sqrt(N_bkg))**2 + sigma_read**2)
  #Hmasked  = np.ma.masked_where((obj.image-N_bkg)<maskval,obj.image-N_bkg)
  #Hmasked2 = np.ma.masked_where((obj.image-N_bkg)<maskval,obj.image)
  #estimate = np.sum(Hmasked)
  #uncertainty = np.sqrt(np.sum(Hmasked2))
  if(display >= 3):
	plt.figure(figsize=(8,11))
        plt.subplot(321)
	plt.imshow(obj.image, interpolation='none', vmin=np.min(obj.image), vmax=np.max(obj.image))
	plt.colorbar()
        plt.title('data')
        plt.subplot(322)
	plt.imshow(obj.moffat_fit_image, interpolation='none', vmin=np.min(obj.image), vmax=np.max(obj.image))
	plt.colorbar()
        plt.title('Moffat Fit')
        plt.subplot(323)
	res = obj.image-obj.moffat_fit_image
	plt.imshow(res, interpolation='none', vmin=-np.max(np.max(res)),  vmax=np.max(np.max(res)))
	plt.colorbar()
	plt.title('residuals')
        plt.subplot(324)
	chi = res/(np.sqrt(obj.image + obj.FM.readnoise**2))
	plt.imshow(chi, interpolation='none', vmin=-np.max(np.abs(chi)), vmax=np.max(np.abs(chi)))
	mx = np.max(np.abs(chi))
	mn = -mx
	plt.colorbar()
	plt.title('chi values')
        plt.subplot(313)
	plt.hist(np.ravel(chi), bins=30)
	plt.xlim(mn,mx)
	plt.suptitle(obj.FM.fits_file_name.split('/')[-1]+'\n'+out, fontsize=20)
	plt.xlabel('chi values')
	plt.subplots_adjust(bottom=0.05)
	plt.title('chi distribution')
        plt.savefig(out+'.png', dpi=50)
  #if(display or uncertainty/estimate>3e-2):
  '''
  if(display>=3):
	  print uncertainty/estimate
	  figure()
	  subplot(221)
          m1 = np.min(obj.image)
          m2 = np.max(obj.image)
          m3 = np.min(obj.image-N_bkg)
          m4 = np.max(obj.image-N_bkg)
	  imshow(obj.image, interpolation='none', vmin=m1, vmax=m2)
	  colorbar()
	  title('Image')
	  subplot(222)
	  imshow(obj.image-N_bkg, interpolation='none', vmin=m3, vmax=m4)
	  colorbar()
	  title('Bkg Subtracted')
	  subplot(223)
	  imshow(Hmasked2, interpolation='none', vmin=m1, vmax=m2)
	  colorbar()
	  subplot(224)
	  imshow(Hmasked, interpolation='none', vmin=m3, vmax=m4)
	  colorbar()
	  show()
  '''
  return estimate, uncertainty


def estimate_total_light_elliptical(obj, N_bkg, sigma_read, display=0, out='out'):
  obj.fit_elliptical_moffat()
  #print obj.moffat_parms
  alpha = obj.elliptical_moffat_parms[1]
  beta = obj.elliptical_moffat_parms[2]
  #print 'a,b',alpha,beta
  # DO SOME COOKIE CUT OUT STUFF LATER, FOR NOW ESTIMATE BY THE NUMBER OF COUNTS DUE TO THE SIGNAL
  estimate = np.sum(obj.elliptical_moffat_fit_image - obj.elliptical_moffat_parms[5]) # magnitude based on value of fitted counts
  uncertainty = np.sqrt(np.sum(obj.elliptical_moffat_fit_image)+sigma_read**2)
  # percent uncertainty is on the level of light.
  #print '\t',N_bkg, obj.moffat_parms[5]
  #maskval = np.sqrt((3.*np.sqrt(N_bkg))**2 + sigma_read**2)
  #maskval = np.sqrt((5.*np.sqrt(obj.moffat_parms[5]))**2 + sigma_read**2)
  #Hmasked  = np.ma.masked_where((obj.image-obj.moffat_parms[5])<maskval,obj.image-obj.moffat_parms[5])
  #Hmasked2 = np.ma.masked_where((obj.image-obj.moffat_parms[5])<maskval,obj.image)
  #estimate = np.sum(Hmasked)
  #frac_uncertainty = np.sqrt(np.sum(Hmasked2)+sigma_read**2) / np.sum(Hmasked)
  print '\timage peak',np.max(obj.image)
  print '\timage bkg',np.min(obj.image)
  #print '\tmasked image sum counts',np.sum(Hmasked2)
  #print '\tmasked bkg subtractected image sum counts',np.sum(Hmasked)
  #print '\tfrac_uncertainty',frac_uncertainty
  print '\tmoffat_integrated count',estimate
  print '\tmoffat estimated uncertainty',uncertainty
  print '\tmoffat estimated frac uncertainty',uncertainty/estimate
  
  #maskval = np.sqrt((3.*np.sqrt(N_bkg))**2 + sigma_read**2)
  #Hmasked  = np.ma.masked_where((obj.image-N_bkg)<maskval,obj.image-N_bkg)
  #Hmasked2 = np.ma.masked_where((obj.image-N_bkg)<maskval,obj.image)
  #estimate = np.sum(Hmasked)
  #uncertainty = np.sqrt(np.sum(Hmasked2))
  if(display >= 3):
	plt.figure(figsize=(8,11))
        plt.subplot(321)
	plt.imshow(obj.image, interpolation='none', vmin=np.min(obj.image), vmax=np.max(obj.image))
	plt.colorbar()
        plt.title('data')
        plt.subplot(322)
	plt.imshow(obj.elliptical_moffat_fit_image, interpolation='none', vmin=np.min(obj.image), vmax=np.max(obj.image))
	plt.colorbar()
        plt.title('Moffat Fit')
        plt.subplot(323)
	res = obj.image-obj.elliptical_moffat_fit_image
	plt.imshow(res, interpolation='none', vmin=-np.max(np.max(res)),  vmax=np.max(np.max(res)))
	plt.colorbar()
	plt.title('residuals')
        plt.subplot(324)
	chi = res/(np.sqrt(obj.image + obj.FM.readnoise**2))
	plt.imshow(chi, interpolation='none', vmin=-np.max(np.abs(chi)), vmax=np.max(np.abs(chi)))
	mx = np.max(np.abs(chi))
	mn = -mx
	plt.colorbar()
	plt.title('chi values')
        plt.subplot(313)
	plt.hist(np.ravel(chi), bins=30)
	plt.xlim(mn,mx)
	plt.suptitle(obj.FM.fits_file_name.split('/')[-1]+'\n'+out, fontsize=20)
	plt.xlabel('chi values')
	plt.subplots_adjust(bottom=0.05)
	plt.title('chi distribution')
        plt.savefig(out+'.png', dpi=50)
  #if(display or uncertainty/estimate>3e-2):
  '''
  if(display>=3):
	  print uncertainty/estimate
	  figure()
	  subplot(221)
          m1 = np.min(obj.image)
          m2 = np.max(obj.image)
          m3 = np.min(obj.image-N_bkg)
          m4 = np.max(obj.image-N_bkg)
	  imshow(obj.image, interpolation='none', vmin=m1, vmax=m2)
	  colorbar()
	  title('Image')
	  subplot(222)
	  imshow(obj.image-N_bkg, interpolation='none', vmin=m3, vmax=m4)
	  colorbar()
	  title('Bkg Subtracted')
	  subplot(223)
	  imshow(Hmasked2, interpolation='none', vmin=m1, vmax=m2)
	  colorbar()
	  subplot(224)
	  imshow(Hmasked, interpolation='none', vmin=m3, vmax=m4)
	  colorbar()
	  show()
  '''
  return estimate, uncertainty

def magnitude(value):
  return -2.5*np.log10(value)

def APASS_zero_points(FM, APASS_table, APASS_rejects, sigma_read, display=0, out = 'out'):
  # KEEP TRACK OF WHICH BAND THE IMAGES ARE AND USE THE CORRESPONDING APASS REFERENCE VALUES.
  filt = 'Sloan_r'
  filt_err = 'r_err'
  if( FM.hdulist[0].header['FILTER'] == 'gp'):
    filt = 'Sloan_g'
    filt_err = 'gerr'
  # MAKE LISTST FOR THE CALIBRATED APASS MAGNITUDES AND THE INSTRUMENT MAGNITUDES
  m_APASS = []
  m_I     = []
  m_APASS_unc = []
  m_I_unc     = []
  alpha = []
  beta  = []
  bkg=[]
  #for k in range(4,6):
  for k in range(0,len(APASS_table)):
    if(APASS_table[filt][k]!='NA' and APASS_table[filt_err][k]!='0' and k not in APASS_rejects): 
	print 'APASS', k
	print 'magnitude %1.2f'%(float(APASS_table[filt][k])) 
	obj = SourceImage(FM, APASS_table['radeg'][k], APASS_table['decdeg'][k], 31)
	N_bkg = np.median(obj.image)
	outTag = out+'_src_%d'%k
	#plt.show()
	try:
		intg, intg_unc = estimate_total_light(obj, N_bkg, sigma_read, display=display, out=outTag)
		#intg, intg_unc = estimate_total_light_elliptical(obj, N_bkg, sigma_read, display=display, out=outTag)

	except:
		continue
	if(obj.fit_ok):
	  m_I.append(magnitude(intg))
	  m_I_unc.append(2.5/2.3*( (intg_unc/intg) - 0.5*(intg_unc/intg)**2 + 1./3.*(intg_unc/intg)**3 ) )
	  m_APASS.append(float(APASS_table[filt][k]))
	  m_APASS_unc.append(float(APASS_table[filt_err][k]))
	  alpha.append(obj.moffat_parms[1])
	  beta.append(obj.moffat_parms[2])
	  bkg.append(obj.moffat_parms[5])
          #print 'k', obj.moffat_parms[1], obj.moffat_parms[2], obj.moffat_parms[5]
	  '''
          figure()
          imshow(obj.image, interpolation='none')
          title('APASS %d, m_APASS %1.2f'%(k,float(APASS_table[filt][k])))
	  '''
  #plt.show()
  # ESTIMATE ZERO POINTS FOR EACH APASS SOURCE
  ZP = np.array(m_APASS)-np.array(m_I)
  # USE ESTIMATED UNCERTINTIES TO WEIGHT THE MEAN AND RMS OF ZERO POINTS
  weight = 1./np.sqrt(np.array(m_I_unc)**2+np.array(m_APASS_unc)**2)
  # ESTIMATE THE MEAN AND WEIGHTED RMS OF THE ZERO POINTS 
  #ZP_mean = sum(ZP)/len(ZP)
  #ZP_rms = np.sqrt(sum(ZP**2)/len(ZP) - ZP_mean**2)
  ZP_mean = sum(weight*ZP)/sum(weight)
  ZP_rms = np.sqrt(sum(weight*(ZP-ZP_mean)**2)/sum(weight))
  print 'Ensable APASS Uncertainties: ZP_rms, m_I_unc, m_APASS_unc',ZP_rms, np.sqrt(np.sum(np.array(m_I_unc)**2/len(m_I))), np.sqrt(np.sum(np.array(m_APASS_unc)**2/len(m_APASS)))
  #print ZP_mean, ZP_rms
  # WEIGHTED RMS
  # ESTIMATE THE WEIGHTED MEAN OF SEEING PARAMETERS
  alpha_mean = sum(weight*alpha)/sum(weight)
  beta_mean  = sum(weight*beta)/sum(weight)

  if(display>0):
  	# SORT BY INCREASING MAGNITUDE FOR EASE OF VIEWING
  	m_APASS_ord,m_I_ord = zip(*sorted(zip(m_APASS, m_I)))
  	m_APASS_ord,m_APASS_unc_ord = zip(*sorted(zip(m_APASS, m_APASS_unc)))
  	m_APASS_ord,m_I_unc_ord = zip(*sorted(zip(m_APASS, m_I_unc)))
  	m_APASS_ord,alpha_ord = zip(*sorted(zip(m_APASS, alpha)))
  	m_APASS_ord,beta_ord = zip(*sorted(zip(m_APASS, beta)))
  	m_APASS_ord,bkg_ord = zip(*sorted(zip(m_APASS, bkg)))
	col='r'
	plt.figure(figsize=(8,12))
	if(filt=='Sloan_r'):
		plt.subplot(511)
		plt.errorbar(m_APASS_ord,np.array(m_APASS_ord)-np.array(m_I_ord),xerr=m_APASS_unc_ord, yerr=np.sqrt(np.array(m_I_unc_ord)**2+np.array(m_APASS_unc_ord)**2) ,fmt='.-', color=col)
		plt.title('Red Filter ZP')
 		plt.ylabel('Zero Point')
		#plt.plot()
	if(filt=='Sloan_g'):
		col='g'
		plt.subplot(511)
		plt.errorbar(m_APASS_ord,np.array(m_APASS_ord)-np.array(m_I_ord), xerr=m_APASS_unc_ord, yerr=np.sqrt(np.array(m_I_unc_ord)**2+np.array(m_APASS_unc_ord)**2),fmt='.-', color=col)
 		plt.ylabel('Zero Point')
 		plt.title('Green Filter ZP')
	plt.xlim(np.min(m_APASS_ord)-0.5, np.max(m_APASS_ord)+0.5)
	plt.grid(True)
	plt.subplot(512)
	plt.plot(m_APASS_ord, alpha_ord, '.-', color=col)
	plt.ylabel('alpha')
	plt.title('Moffat alpha')
	plt.grid(True)
	plt.xlim(np.min(m_APASS_ord)-0.5, np.max(m_APASS_ord)+0.5)
	plt.subplot(513)
	plt.plot(m_APASS_ord, beta_ord, '.-', color=col)
	plt.ylabel('beta')
	plt.title('Moffat beta')
	plt.grid(True)
	plt.xlim(np.min(m_APASS_ord)-0.5, np.max(m_APASS_ord)+0.5)
	plt.subplot(514)
	plt.plot(m_APASS_ord, moffat_fwhm(np.array(alpha_ord), np.array(beta_ord)), '.-', color=col)
	plt.ylabel('FWHM')
	plt.title('Moffat FWHM')
	plt.grid(True)
	plt.xlim(np.min(m_APASS_ord)-0.5, np.max(m_APASS_ord)+0.5)
	plt.subplot(515)
	plt.plot(m_APASS_ord, bkg_ord, '.-', color=col)
	plt.title('Fitted Nbkg')
	plt.ylabel('Counts')
	plt.grid(True)
	plt.subplots_adjust(left=0.2, right=0.95, hspace=0.4)
	plt.xlabel('APASS Source Magnitude', fontsize=14)
	plt.xlim(np.min(m_APASS_ord)-0.5, np.max(m_APASS_ord)+0.5)
	plt.suptitle(FM.fits_file_name.split('/')[-1], fontsize=18)
	plt.savefig(out+'.png', dpi=50)
  return ZP_mean, ZP_rms, alpha_mean, beta_mean


def quadFit(FM, ra_qsr, dec_qsr, ZP_mean, ZP_rms, alpha, beta, N_px, outputFileTag='out'):
  # GET THE QUASAR IMAGE
  obj = SourceImage(FM, ra_qsr, dec_qsr, N_px)

  N_bkg = np.median(obj.image)
  intg, intg_unc = estimate_total_light(obj,N_bkg, FM.readnoise, display=False)

  # DETERMINE THE LEFT RIGHT ORIENTATION OF THE IMAGE BASED ON ITS RIGHT ASCENCION DIFFERENCE WITH PIXEL VALUES
  xv,yv = FM.bigw.wcs_world2pix(ra_qsr,dec_qsr,1)
  r0,d0 = FM.bigw.wcs_pix2world(xv,yv,1)
  r1,d1 = FM.bigw.wcs_pix2world(xv+1,yv+1,1)
  fl = True
  if(r1-r0<0): fl =False

  # SET AN INITIAL AMPLITUDE SCALE, THESE VALUES ARE APPROXIMATED FROM Kochenek, 2005 
  amp_scale = np.max(obj.image)
  x0,y0 = 0., 0.
  amp0 = amp_scale
  amp1 = amp_scale/1.2
  amp2 = amp_scale/1.5
  amp3 = amp_scale/2.
  # DEFINE ARRAY OF INPUT PARAMETERS
  # PRODUCE A MODELED IMAGE
  qim = obj.quad_image_model(x0,y0,amp0,amp1,amp2,amp3, alpha, beta, N_bkg, len(obj.image), fl)
  # CROSS CORRELATE THE DATA TO THE MODEL TO FIND OUT WHERE IT IS ON THE PIXEL GRID
  corr = signal.correlate2d(obj.image, qim, boundary='symm', mode='same')
  # GET THE LOCATION OF THE CORRELATION PEAK
  corr_x_max, corr_y_max = np.unravel_index(corr.argmax(), corr.shape)
  # GET THE IMAGE LOCATION CORRECTIONS
  dx = corr_x_max-15
  dy = corr_y_max-15
  # SOME FUNINESS HERE. I DETERMINED IT FROM SCANNING THROUGH SEVERAL IMAGES
  if(corr_x_max<15 and corr_y_max>=15):
    dy *= -1
  # DEFINE NEW PARAMETERS FOR THE MODELED IMAGE
  x0 += dx
  y0 += dy
  amp0*=np.max(obj.image)/np.max(qim)
  amp1*=np.max(obj.image)/np.max(qim)
  amp2*=np.max(obj.image)/np.max(qim)
  amp3*=np.max(obj.image)/np.max(qim)
  # PRODUCE THE MODEL IMAGE
  qim2 = obj.quad_image_model(x0,y0,amp0,amp1,amp2,amp3, alpha, beta, N_bkg, 31, fl)

  # PRODUCE THE PARAMETERS FOR THE MODEL IMAGE
  theta2 = [x0,y0,amp0,amp1,amp2,amp3, alpha, beta, N_bkg]
  #theta2 = [x0,y0,amp0,amp1,amp2,amp3, N_bkg]

  # ESTIMATE ITS DISTRIBUTION OF CHI VALUES
  chi2 = obj.moffat_chi_vals(theta2, 31, fl)
  #chi2 = obj.moffat_chi_vals2(theta2, alpha, beta,31, fl)

  # FIT A MODEL TO THE DATA
  results = opt.minimize(obj.moffat_chi_sq, theta2,  args=(31, fl), method='Nelder-Mead')
  #results = opt.minimize(obj.moffat_chi_sq2, theta2,  args=(alpha, beta, 31, fl), method='Nelder-Mead')

  # I FIND IT PAYS TO KEEP THE MINIMIZATION GOING FOR A FEW MORE ROUNDS
  for k in range(0,4):
	results = opt.minimize(obj.moffat_chi_sq, results.x,  args=(31, fl), method='Nelder-Mead')
	#results = opt.minimize(obj.moffat_chi_sq2, results.x,  args=(alpha, beta, 31, fl), method='Nelder-Mead')
  # GET THE FITTED IMAGE MODEL
  x0,y0,amp0,amp1,amp2,amp3, alpha, beta, N_bkg = results.x
  #x0,y0,amp0,amp1,amp2,amp3, N_bkg = results.x
  qim3 = obj.quad_image_model(x0,y0,amp0,amp1,amp2,amp3, alpha, beta, N_bkg, 31, fl)

  # ESTIMATE ITS CHI VALUES, ITS PEAK CHI VALUE AND LOCATION ON THE MAP
  chi_vals = obj.moffat_chi_vals(results.x, 31, fl)
  #chi_vals = obj.moffat_chi_vals2(results.x, alpha, beta,31, fl)

  chi_x_max, chi_y_max = np.unravel_index(np.abs(chi_vals).argmax(), chi_vals.shape)
  #print 'alpha, beta, N_bkg', alpha, beta, N_bkg
  #print 'max chi x,y, val', chi_x_max, chi_y_max, np.max(np.abs(chi_vals))
  #print results.x

  xg, yg = np.mgrid[:31, :31]

  parms1 = [results.x[2], results.x[6], results.x[7], 15, 15, results.x[8]]
  moffat_fit_image_1  = twoD_Moffat((xg, yg), *parms1).reshape(31,31)

  parms2 = [results.x[3], results.x[6], results.x[7], 15, 15, results.x[8]]
  moffat_fit_image_2  = twoD_Moffat((xg, yg), *parms2).reshape(31,31)

  parms3 = [results.x[4], results.x[6], results.x[7], 15, 15, results.x[8]]
  moffat_fit_image_3  = twoD_Moffat((xg, yg), *parms3).reshape(31,31)

  parms4 = [results.x[5], results.x[6], results.x[7], 15, 15, results.x[8]]
  moffat_fit_image_4  = twoD_Moffat((xg, yg), *parms4).reshape(31,31)

  print '\t',parms1
  print '\t',parms2
  print '\t',parms3
  print '\t',parms4
  
  FM.qsr_min_parms = results.x
  '''
  plt.figure()
  plt.subplot(221)
  plt.imshow(moffat_fit_image_1-results.x[8])
  plt.title(np.sum(moffat_fit_image_1-results.x[8]))
  plt.colorbar()
  plt.subplot(222)
  plt.imshow(moffat_fit_image_2-results.x[8])
  plt.title(np.sum(moffat_fit_image_2-results.x[8]))
  plt.colorbar()
  plt.subplot(223)
  plt.imshow(moffat_fit_image_3-results.x[8])
  plt.title(np.sum(moffat_fit_image_3-results.x[8]))
  plt.colorbar()
  plt.subplot(224)
  plt.imshow(moffat_fit_image_4-results.x[8])
  plt.title(np.sum(moffat_fit_image_4-results.x[8]))
  plt.colorbar()
  plt.show()
  '''

  '''
  a1 = moffat_analytical_integral(results.x[2], alpha, beta)
  a2 = moffat_analytical_integral(results.x[3], alpha, beta)
  a3 = moffat_analytical_integral(results.x[4], alpha, beta)
  a4 = moffat_analytical_integral(results.x[5], alpha, beta)
  '''
  print '\tpeaks', np.sum(moffat_fit_image_1-results.x[8]), np.sum(moffat_fit_image_2-results.x[8]), np.sum(moffat_fit_image_3-results.x[8]), np.sum(moffat_fit_image_4-results.x[8])
  a1 = np.sum(moffat_fit_image_1-results.x[8])
  a2 = np.sum(moffat_fit_image_2-results.x[8])
  a3 = np.sum(moffat_fit_image_3-results.x[8])
  a4 = np.sum(moffat_fit_image_4-results.x[8])

  e1 = np.sqrt( np.sum(moffat_fit_image_1) + FM.readnoise**2 )
  e2 = np.sqrt( np.sum(moffat_fit_image_2) + FM.readnoise**2 )
  e3 = np.sqrt( np.sum(moffat_fit_image_3) + FM.readnoise**2 )
  e4 = np.sqrt( np.sum(moffat_fit_image_4) + FM.readnoise**2 )

  print '\ta1', a1
  print '\te1', e1
  print '\ta2', a2
  print '\te2', e2
  print '\ta3', a3
  print '\te3', e3
  print '\ta4', a4
  print '\te4', e4
  print ''

  m1 = magnitude(a1) + ZP_mean
  m2 = magnitude(a2) + ZP_mean
  m3 = magnitude(a3) + ZP_mean
  m4 = magnitude(a4) + ZP_mean

  me1 = 2.5/2.3*(e1/a1-1./2.*(e1/a1)**2 + 1./3.*(e1/a1)**3)
  me2 = 2.5/2.3*(e2/a2-1./2.*(e2/a2)**2 + 1./3.*(e2/a2)**3)
  me3 = 2.5/2.3*(e3/a3-1./2.*(e3/a3)**2 + 1./3.*(e3/a3)**3)
  me4 = 2.5/2.3*(e4/a4-1./2.*(e4/a4)**2 + 1./3.*(e4/a4)**3)

  print '\tm1', m1
  print '\tme1', me1
  print '\tm2', m2
  print '\tme2', me2
  print '\tm3', m3
  print '\tme3', me3
  print '\tm4', m4
  print '\tme4', me4

  me1 = np.sqrt(me1**2 + ZP_rms**2 )
  me2 = np.sqrt(me2**2 + ZP_rms**2 )
  me3 = np.sqrt(me3**2 + ZP_rms**2 )
  me4 = np.sqrt(me4**2 + ZP_rms**2 )

  print ''
  print '\tm1', m1
  print '\tme1', me1
  print '\tm2', m2
  print '\tme2', me2
  print '\tm3', m3
  print '\tme3', me3
  print '\tm4', m4
  print '\tme4', me4

  #print m1,me1, m2,me2, m3, me3, m4, me4

  
  plt.figure(figsize=(12, 12)) 
  plt.subplot(331)
  plt.imshow(obj.image, interpolation='none')
  plt.colorbar()
  plt.contour(obj.image, colors='k')
  plt.title('data')
  plt.subplot(332)
  plt.imshow(qim, interpolation='none')
  plt.colorbar()
  plt.contour(obj.image, colors='k')
  plt.title('first stab\nw/ data contour')
  plt.subplot(333)
  plt.imshow(corr, interpolation='none')
  plt.colorbar()
  plt.title('first stab & data\ncross-correlation')
  plt.subplot(334)
  plt.imshow(qim2, interpolation='none')
  plt.colorbar()
  plt.contour(obj.image, colors='k')
  plt.title('initial model\nw/ data controur')
  plt.subplot(335)
  mv = np.max(np.abs(obj.image - qim2))
  plt.imshow(obj.image - qim2, cmap='jet', interpolation='none', vmin=-mv, vmax=mv)
  plt.colorbar()
  plt.contour(obj.image, colors='k')
  plt.title('initial model residuals')
  plt.subplot(336)
  mv = np.max(np.abs(chi2))
  plt.imshow(chi2, cmap='jet', interpolation='none', vmin=-mv, vmax=mv)
  plt.colorbar()
  plt.contour(obj.image, colors='k')
  plt.subplot(337)
  plt.imshow(qim3, interpolation='none', vmin=np.min(obj.image), vmax=np.max(obj.image))
  plt.colorbar()
  plt.contour(obj.image, colors='k')
  plt.title('fitted model \n w/ data contour')
  plt.subplot(338)
  mv = np.max(np.abs(obj.image - qim3))
  plt.imshow(obj.image - qim3, cmap='jet', interpolation='none', vmin=-mv, vmax=mv)
  plt.colorbar()
  plt.contour(obj.image, colors='k')
  plt.title('fitted model residuals')
  plt.subplot(339)
  mv = np.max(np.abs(chi_vals))
  plt.imshow(chi_vals, cmap='jet', interpolation='none', vmin=-mv, vmax=mv)
  plt.colorbar()
  plt.contour(obj.image, colors='k')
  plt.plot([chi_y_max], [chi_x_max],'o', mfc='none', mec='k', mew=3, ms=15)
  plt.title('fitted Model Chi vals')
  plt.suptitle(obj.FM.fits_file_name.split('/')[-1],fontsize=20)
  plt.savefig(outputFileTag+'_quadFit.png', dpi=50)

  plt.figure()
  plt.hist(chi_vals.flat)
  plt.xlabel('chi values')
  plt.ylabel('counts')
  plt.suptitle(obj.FM.fits_file_name.split('/')[-1],fontsize=20)
  plt.savefig(outputFileTag+'_quadFitChi.png', dpi=50)
  #plt.savefig()

  flat_chiVals = np.ravel(chi_vals)
  chiSq= np.sum(flat_chiVals**2)/(len(flat_chiVals)-len(results.x))
  #chiSq = np.sum(chi_vals.flat**2)/(len(chi_vals.flat)-len(results.x))
  
  #print chiSq
  #exit()
  return m1, me1, m2, me2, m3, me3, m4, me4, chiSq, np.max(np.abs(chi_vals))
  #plt.show()
  #figure()
  #hist(chi_vals)
  #show()

  '''
  if(count2==11):
	plt.figure()
	estimate_total_light(obj, N_bkg, sigma_read, display=True)
	plt.figure()
	levels=arange(0,np.max(obj.image),100)[::-1]
	plt.imshow(obj.image-N_bkg, interpolation='none')
	plt.colorbar()
	plt.show()
  '''

#def emceeQuadFit(FM, ra_qsr, dec_qsr, ZP_mean, ZP_rms, alpha, beta, m1, m2, m3, m4, me4, N_px, outputFileTag='out'):
def emceeQuadFit(FM, ra_qsr, dec_qsr, ZP_mean, ZP_rms, alpha, beta, m1, m2, m3, m4, N_px, outputFileTag='out'):
  print '\n\t####################################'
  print   '\t# STARTING EMCEE ###################'
  print   '\t####################################\n'
  obj = SourceImage(FM, ra_qsr, dec_qsr, N_px)
  N_bkg_guess = np.median(obj.image)

  x0,y0, amp1, amp2, amp3, amp4, alpha_qsr, beta_qsr, N_bkg_qsr =FM.qsr_min_parms
  #amp1 = 10**((m1-ZP_mean)/(-2.5))
  #amp2 = 10**((m2-ZP_mean)/(-2.5))
  #amp3 = 10**((m3-ZP_mean)/(-2.5))
  #amp4 = 10**((m4-ZP_mean)/(-2.5))
  #N_bkg = np.median(obj.image)

  #obj = SourceImage(FM, ra_qsr, dec_qsr, N_px)

  def lnprior(_theta):
    _x0, _y0, \
    _amplitude0, \
    _amplitude1, \
    _amplitude2, \
    _amplitude3, \
    _alpha, \
    _beta, \
    _N_bkg = _theta
    #print 'N_bkg_guess', N_bkg_guess
    #exit()
    if  np.abs(_x0)<N_px \
	and np.abs(_y0)<N_px \
	and _amplitude0>0. and _amplitude0<100.*amp1\
	and _amplitude1>0. and _amplitude1<100.*amp2\
	and _amplitude2>0. and _amplitude2<100.*amp3\
	and _amplitude3>0. and _amplitude3<100.*amp4\
	and _alpha>0. and _alpha < 100.*alpha\
	and _beta>0. and  _beta  < 100.*beta \
	and _N_bkg>0. and _N_bkg < 100.*N_bkg_guess:
	  return 0.0  
    else:
        return -np.inf

  def lnprob(theta, obj, N_pix, flip):
    lp = lnprior(theta)
    #print lp, theta
    if not np.isfinite(lp):
        return -np.inf
    mll = obj.moffat_chi_sq(theta, N_pix, flip)
    #print mll
    if not np.isfinite(mll):
        return -np.inf
    #print '%1.2e\t%1.2e'%(lp,mll), theta
    return lp - mll

  # DETERMINE THE IMAGE ORIENTATION
  xv,yv = FM.bigw.wcs_world2pix(ra_qsr,dec_qsr,1)
  r0,d0 = FM.bigw.wcs_pix2world(xv,yv,1)
  r1,d1 = FM.bigw.wcs_pix2world(xv+1,yv+1,1)
  fl = True
  if(r1-r0<0): fl =False

  #################################################
  #################################################
  #################################################

  # RUN MC!
  theta = FM.qsr_min_parms
  #theta = [x0,y0,amp1,amp2,amp3,amp4, alpha, beta, N_bkg]

  print theta

  ndim =len(theta)
  print 'ndim', ndim
  nwalkers = 100
  n_burn_in_iterations = 10
  n_iterations = 10000

  
  prior_vals=theta
  pos = [prior_vals + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
  r=np.random.randn(ndim)
  sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(obj,N_px, fl))

  #print 'burn-in sampler.run_mcmc'
  sampler.run_mcmc(pos, n_burn_in_iterations)

  samples = sampler.chain[:, int(0.5*float(n_burn_in_iterations)):, :].reshape((-1, ndim))
  print("\tMean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
  parms = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))
  #labels = [ r"$x_0$", r"$y_0$", "amp1", "amp2", "amp3", "amp4", r"$\alpha$", "$\beta$", "N$_{bkg}$"]
  labels = [ "x0", "y0", "amp1", "amp2", "amp3", "amp4", "alpha", "beta", "N_bkg"]
  fig= triangle.corner(samples, labels=labels)
  count=0
  for ax in fig.get_axes():
    count+=1
    #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
    #ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
    if((count-1)%ndim==0 ): 
	ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
	ax.xaxis.set_label_coords( 0.5, -0.5)
	ax.yaxis.set_label_coords(-0.5, 0.5)

    if(count>=ndim*(ndim-1)): 
	ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
	ax.xaxis.set_label_coords( 0.5, -0.5)
	ax.yaxis.set_label_coords(-0.5,  0.5)

  plt.subplots_adjust(bottom=0.075, left=0.075)
  fig.savefig("%s_burn_in_triangle.png"%(outputFileTag))

  for k in range(0,len(parms)):
    print '%s\t%+1.2e\t%+1.2e\t%+1.2e\t%+1.2e'%(labels[k], theta[k], parms[k][0],parms[k][1],parms[k][2])

  #exit()
  
  sampler.reset()
  sampler.run_mcmc(pos, n_iterations)
  samples = sampler.chain[:, int(0.5*float(n_iterations)):, :].reshape((-1, ndim))
  print("\tMean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
  parms = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))

  #print parms.shape
  #print parms[k][0]

  for k in range(0,len(parms)):
    print '%s\t%+1.2e\t%+1.2e\t%+1.2e\t%+1.2e'%(labels[k], theta[k], parms[k][0],parms[k][1],parms[k][2])

  fig= triangle.corner(samples, labels=labels)
  count=0
  for ax in fig.get_axes():
    count+=1
    #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
    #ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
    if((count-1)%ndim==0 ): 
	ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
	ax.xaxis.set_label_coords( 0.5, -0.5)
	ax.yaxis.set_label_coords(-0.5, 0.5)

    if(count>=ndim*(ndim-1)): 
	ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
	ax.xaxis.set_label_coords( 0.5, -0.5)
	ax.yaxis.set_label_coords(-0.5,  0.5)

  plt.subplots_adjust(bottom=0.075, left=0.075)
  fig.savefig("%s_triangle.png"%(outputFileTag))

  np.savez(outputFileTag+'_chains.npz', sampler.chain[:, 0:, :].reshape((-1, ndim)))
  print 'chain saving complete'




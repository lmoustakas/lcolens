import numpy as np
import datetime
import os
from astropy.io import fits
from astropy.io import ascii
from astropy import wcs
#from astroscrappy import detect_cosmics

#from wcsaxes import WCSAxes
#import aplpy
import pylab as plt
import scipy.optimize as opt
from scipy.optimize import minimize
from scipy.optimize import leastsq
from scipy.special import genlaguerre
from scipy.misc import factorial
import matplotlib
from scipy import signal
import emcee
import corner
import psf
#from astroscrappy import detect_cosmics
import integrate
integrate = integrate.integrate()
import time
from matplotlib.colors import LogNorm

################################################################

def HMS2deg(ra='', delimiter=':'):
    RA, rs = '', 1  
    if ra:
        H, M, S = [float(i) for i in ra.split(delimiter)]
        if str(H)[0] == '-':
            rs, H = -1, abs(H)
        deg = (H*15) + (M/4) + (S/240)
        RA = '{0}'.format(deg*rs)
    return RA

################################################################

def DMS2deg(dec='', delimiter=':'):
    DEC, ds = '', 1
    if dec:
        D, M, S = [float(i) for i in dec.split(delimiter)]
        if str(D)[0] == '-':
            ds, D = -1, abs(D)
        deg = D + (M/60) + (S/3600)
        DEC = '{0}'.format(deg*ds)
    return DEC


################################################################

def convertRaAndDec(ra, dec):
    if (-1 != ra.find(':')):
        ra = HMS2deg(ra, ':')
    elif (-1 != ra.find(' ')):
        ra = HMS2deg(ra, ' ')
    if (-1 != dec.find(':')):
        dec = DMS2deg(dec, ':')
    elif (-1 != dec.find(' ')):
        dec = DMS2deg(dec, ' ')
    return ra, dec

################################################################

def time_order_fits_files(data_dir):
	# GET ALL DATA FILE NAMES IN DATA DIRECTORY AND SORT THEM BY MJD
	print 'READING THE FITS FILE OBSERVATION TIMES FOR ORDERING PURPOSES'
	fnms = os.listdir(data_dir)
	print 'Number of files in directory:', len(fnms)
	fnms_filtered = []
	for f in fnms:
	  #print f[-5:], f
	  if(len(f.split('-'))>3 and f[-5:]=='.fits'):
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

    #print "self.hdulist.info()"
    #print self.hdulist.info()
    '''
    print "self.hdulist[0].header"
    print self.hdulist[0].header
    print "self.hdulist[1].header"
    print self.hdulist[1].header
    print "self.hdulist[2].header"
    print self.hdulist[2].header
    print "self.hdulist[3].header"
    print self.hdulist[3].header
    '''
    #exit() 
    #self.hdulist.info()
    # this slows down the code significantly!!!!
    # make sure masking is an option whose default is false.
    #mask, clean = detect_cosmics(self.hdulist[0].data, sigfrac=0.15, sigclip=4, objlim=4, cleantype='idw')
    #self.hdulist[0].data = np.ma.masked_where(mask==1,self.hdulist[0].data)
    #self.hdulist[0].data = np.ma.masked_where(self.hdulist[0].data<0.,self.hdulist[0].data)
    self.image_data = self.hdulist[0].data

    # MULTIPLY BY THE GAIN
    # this has been commented out. We model the noise variance for CCD values throughout. 
    # self.image_data *= float(self.hdulist[0].header['GAIN'])
    self.bigw=wcs.WCS(self.hdulist[0].header)

    print '\nFITSmanager'
    print '\tGAIN'.ljust(9), '=', self.hdulist[0].header['GAIN']
    print '\tPIXSCALE'.ljust(9), '=', self.hdulist[0].header['PIXSCALE']

  def estimate_read_noise(self, display=0, out = 'out'):
	print '\nestimate_read_noise'
	readnoise=[]
	gain = self.hdulist[0].header['GAIN']
	for k in range(1000):
		#print self.image_data.shape
		x = np.random.randint(100, self.image_data.shape[0]-100)
		y = np.random.randint(100, self.image_data.shape[1]-100)
		r, d = self.bigw.wcs_pix2world(x,y, 1)
		#print 'x,y',x,y
		#print 'ra,dec',r,d
		img = self.sub_stamp(r, d, 31)
		#img/=float(self.hdulist[0].header['GAIN']) # remove gain correction. read noise is on electrons, not photons ?		
		count = 0
		while(np.max(img)>np.median(img)+5.*np.sqrt(np.median(img)) or np.min(img)<np.median(img)-5.*np.sqrt(np.median(img))):
			count+=1
	  		#print 'trial', count
			x = np.random.randint(100, self.image_data.shape[0]-100)
			y = np.random.randint(100, self.image_data.shape[1]-100)
			r, d = self.bigw.wcs_pix2world(x,y, 1)
			#print 'x,y',x,y
			#print 'ra,dec',r,d
			img = self.sub_stamp(r, d, 31)
		
		h,b = np.histogram(img.flat, bins=30)
		x=b[:-1]+(b[1]-b[0])/2.
		try:
			popt, pcov = opt.curve_fit(gauss_1d, (x), h, p0=[np.max(h), np.sqrt(np.median(x)), np.median(x)])
		except:
			continue
		if(popt[0]<=0. or popt[2]<0):
			continue
		if(popt[1]**2 - popt[2]/gain <= 0.):
			continue
		rn = np.sqrt(popt[1]**2 - popt[2]/gain)
		if(rn == rn):
			readnoise.append(rn)
	 		#print k,'readnoise', popt[1]-np.sqrt(popt[0])
	 		#print 'median', np.median(img)
			#print popt
	self.readnoise = np.median(readnoise)
	print '\tmedian read noise = %1.2f'%self.readnoise
	#print np.median(readnoise)
	if(display>=2):
		plt.figure(figsize=(7,7))
		plt.subplot(221)
		plt.imshow(img, interpolation='none', cmap='viridis')
		plt.colorbar()
		plt.title('Random Image\nSub-Stamp')
		plt.xlabel('pixel')
		plt.ylabel('pixel')
		plt.subplot(222)
		plt.step(x,h, lw=2)
		plt.plot(x,gauss_1d(x,*popt), 'r--', lw=2) 
		plt.title('Distribution of Counts\nw/ Gaussian Fit')
		plt.xlabel('CCD Counts')
		plt.xticks(rotation=45)
		plt.subplot(212)
		a,b,c = plt.hist(readnoise, bins=int(np.max(readnoise))+5, range=(0,np.max(readnoise)+5))
		plt.plot([np.median(readnoise), np.median(readnoise)],[0., 1.2*max(a)],'r-', lw=2, label='median')
		plt.legend(loc=1)
		plt.ylim(0.,1.2*max(a))
		plt.xlabel('System Noise, CCD counts')
		plt.title('Multiple System Noise Estimates')
		plt.suptitle(self.fits_file_name.split('/')[-1], fontsize=20)
		plt.subplots_adjust(top=0.85, bottom=0.1, wspace=0.3, hspace=0.45, left=0.1, right=0.95)
		plt.savefig(out+'.png', dpi=50)
		#plt.savefig(out+'.pdf')
	#return self.readnoise
	
################################################################

  def histogram_image_values(self,NBINS=5000, rng_max=200000.):
    self.hist, self.bin_edges = np.histogram(self.image_data.flat, bins=NBINS, range=(0.,rng_max))
    self.BW = self.bin_edges[1] - self.bin_edges[0]

################################################################

  def plot_image_values(self,NBINS=5000, rng_max=200000.):
    #histogram = plt.hist(self.image_data.flat, log=True, bins=NBINS)
    self.histogram_image_values(NBINS=NBINS, rng_max=rng_max)
    plt.plot(self.bin_edges[:-1]+self.BW/2.,self.hist+1.e-1, label=self.fits_file_name)
    plt.xlabel('Image Value')
    plt.ylabel('Counts')

################################################################

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

################################################################

  def isFlipped(self, ra, dec):
    xv,yv = self.bigw.wcs_world2pix(ra,dec,1)
    r0,d0 = self.bigw.wcs_pix2world(xv,yv,1)
    r1,d1 = self.bigw.wcs_pix2world(xv+1,yv+1,1)
    #print 'flip', r1-r0, d1-d0
    if(r1-r0>0.): 
	return True
    return False

################################################################

  def plot_image_movie(self, ra, dec, ra_images, dec_images, ax1, ax2, ax3, ZP_flx, theta, dx=0., dy=0.,Npx=31, out='out'):
    # GET THE QUASAR IMAGE
    x0, y0, amp1, amp2, amp3, amp4, alpha, beta, nbkg = theta
    fl = self.isFlipped(ra,dec)
	# convert quasar and quasar image positions to pixel locations
    x,y = self.bigw.wcs_world2pix(ra,dec,1) 
    print ra_images, dec_images
    dec_images = dec + dec_images/3600.
    ra_images = ra   + ra_images/3600./np.cos(dec*np.pi/180.)
    x_images,y_images = self.bigw.wcs_world2pix(ra_images,dec_images,1)
    x_images -= x - (Npx-1)/2
    y_images -= y - (Npx-1)/2
    print x, x_images
    print y, y_images

	# get data and model images
    obj = SourceImage(self, ra, dec, Npx)
    qim = obj.quad_image_model(x0, y0, x_images, y_images, amp1, amp2, amp3, amp4, alpha, beta, nbkg, Npx, fl)
    #qim = obj.quad_image_model(0., 0., x_images, y_images, amp1, amp2, amp3, amp4, alpha, beta, nbkg, Npx, fl)
    '''
    plt.figure()
    plt.subplot(121)
    plt.imshow(obj.image, interpolation='none')
    plt.contour(obj.image, colors='k')
    plt.subplot(122)
    plt.imshow(qim, interpolation='none')
    plt.contour(obj.image)
    #plt.contour(self.image_data[y-(Npx-1)/2:y+(Npx-1)/2+1,x-(Npx-1)/2:x+(Npx-1)/2+1])
    print x,y,x0,y0
    plt.savefig('tmp.png')
	'''


    x_img, y_img = obj.quad_image_pix_ref(x0, y0, x_images, y_images, fl)
    x_img -= (Npx-1)/2+1
    y_img -= (Npx-1)/2+1
    yg, xg =      np.mgrid[y-(Npx-1)/2:y+(Npx-1)/2+1,x-(Npx-1)/2:x+(Npx-1)/2+1]
    r, d = self.bigw.wcs_pix2world(xg-y0,yg-x0, 1)
    d = (d-dec)*3600.
    r = (r-ra) * np.cos(dec*np.pi/180.)*3600.
    fwhm = alpha*2.*np.sqrt(2.**(1/beta)-1.)*float(self.hdulist[0].header['PIXSCALE'])

    #vals = self.image_data[y-(Npx-1)/2:y+(Npx-1)/2+1,x-(Npx-1)/2:x+(Npx-1)/2+1]
    #vals = self.image_data[y-(Npx-1)/2:y+(Npx-1)/2+1,x-(Npx-1)/2:x+(Npx-1)/2+1]
    vals = obj.image
    min_ra = -4.
    max_ra = 4.
    min_dec = -4.
    max_dec = 4.
    vals -= nbkg
    vals*=ZP_flx

    #template = 0.*vals.copy()
    #for i in range(0,4):
    #    template[int(np.rint(y_im[i])-y+(Npx-1)/2),int(np.rint(x_im[i])-x+(Npx-1)/2)]=3.e-9
    vals0 = vals.copy()
    #print xg.shape, yg.shape, vals.shape
    ax1 = plt.subplot(3,2,2)
    #plt.pcolor(xg-x, yg-y, vals, cmap='jet', vmin=0., vmax=3.8e-9)
    plt.pcolor(r, d, vals, cmap='jet', vmin=0., vmax=3.8e-9)
    cb = plt.colorbar()
    plt.plot([min_ra,max_ra], [0.,0.], 'k--')
    plt.plot([0.,0.], [min_dec,max_dec], 'k--')

    plt.axis('equal')
    plt.xlim(max_ra, min_ra)
    plt.ylim(min_dec, max_dec)
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax1.yaxis.set_major_formatter(y_formatter)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax1.xaxis.set_major_formatter(x_formatter)
    #plt.xticks(rotation=30)
    #plt.title('x0 %1.3f y0 %1.3f'%(x0,y0))
    #plt.title('Data')
    plt.xlabel('$\Delta$RA,  arcseconds')
    plt.ylabel('$\Delta$Dec, arcseconds')
    ax1.text( 3.8,  3.5, 'Data', fontsize=30, color='white', fontweight='bold')
    ax1.text( 3.0,  0.0, 'M1', fontsize=30, color='white', fontweight='bold')
    ax1.text(-2.0, -0.5, 'M2', fontsize=30, color='white', fontweight='bold')
    ax1.text( 0.0,  2.0, 'S1', fontsize=30, color='white', fontweight='bold')
    ax1.text(0.5, -2.75, 'S2', fontsize=30, color='white', fontweight='bold')

    ax2 = plt.subplot(3,2,4)
    
    qim -= nbkg
    qim *= ZP_flx
    #plt.pcolor(xg-x,yg-y,qim, cmap='jet', vmin=0., vmax=3.8e-9)
    plt.pcolor(r,d,qim, cmap='jet', vmin=0., vmax=3.8e-9)
    cb = plt.colorbar()
    plt.plot([min_ra,max_ra], [0.,0.], 'k--')
    plt.plot([0.,0.], [min_dec,max_dec], 'k--')
    plt.axis('equal')
    plt.xlim(max_ra, min_ra)
    plt.ylim(min_dec, max_dec)
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax2.yaxis.set_major_formatter(y_formatter)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax2.xaxis.set_major_formatter(x_formatter)
    #plt.xticks(rotation=30)
    plt.xlabel('$\Delta$RA,  arcseconds')
    plt.ylabel('$\Delta$Dec, arcseconds')
    #plt.title('Model')
    ax2.text(3.8, 3.5, 'Model', fontsize=30, color='white', fontweight='bold')

    circle1=plt.Circle((3,-3),fwhm/2.,color='white', fill=False, linewidth=3, linestyle='dashed')
    plt.gcf().gca().add_artist(circle1)
    #plt.circles(3, -3, fwhm/2., alpha=1, lw=5, edgecolor='white', transform=ax2.transAxes)

    ax3 = plt.subplot(3,2,6)
    #plt.pcolor(xg-x,yg-y,vals-qim)
    plt.pcolor(r,d,vals-qim, vmin=-5.5e-10, vmax=5.5e-10)
    cb = plt.colorbar()
    plt.plot([min_ra,max_ra], [0.,0.], 'k--')
    plt.plot([0.,0.], [min_dec,max_dec], 'k--')
    plt.axis('equal')
    plt.xlim(max_ra, min_ra)
    plt.ylim(min_dec, max_dec)
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax3.yaxis.set_major_formatter(y_formatter)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax3.xaxis.set_major_formatter(x_formatter)
    #plt.xticks(rotation=30)
    plt.xlabel('$\Delta$RA,  arcseconds')
    plt.ylabel('$\Delta$Dec, arcseconds')
    ax3.text(3.8, 3.5, 'Residuals', fontsize=30, color='black', fontweight='bold')

    #plt.title('Residuals')

################################################################

  def plot_image_movie_lens_gal(self, ra, dec, ra_images, dec_images, ra_lens, dec_lens, ax1, ax2, ax3, ZP_flx, theta, dx=0., dy=0.,Npx=31, out='out'):
    # GET THE QUASAR IMAGE
    x0, y0, amp1, amp2, amp3, amp4, amp_lens, r0_lens, q_lens, posang_lens,  alpha, beta, nbkg = theta
    fl = self.isFlipped(ra,dec)
        # convert quasar and quasar image positions to pixel locations
    x,y = self.bigw.wcs_world2pix(ra,dec,1)
    print ra_images, dec_images
    dec_images = dec + dec_images/3600.
    ra_images = ra   + ra_images/3600./np.cos(dec*np.pi/180.)
    x_images,y_images = self.bigw.wcs_world2pix(ra_images,dec_images,1)
    x_images -= x - (Npx-1)/2
    y_images -= y - (Npx-1)/2
    x_lens,y_lens = self.bigw.wcs_world2pix(ra_lens,dec_lens,1)
    x_lens -= x - (Npx-1)/2
    y_lens -= y - (Npx-1)/2
    print x, x_images
    print y, y_images

        # get data and model images
    obj = SourceImage(self, ra, dec, Npx)
    qim = obj.quad_image_model(x0, y0, x_images, y_images, amp1, amp2, amp3, amp4, alpha, beta, nbkg, Npx, fl, x_lens=x_lens, y_lens=y_lens, amp_lens=amp_lens, r0_lens=r0_lens, q_lens=q_lens, posang_lens=posang_lens)
    #qim = obj.quad_image_model(0., 0., x_images, y_images, amp1, amp2, amp3, amp4, alpha, beta, nbkg, Npx, fl)
    '''
    plt.figure()
    plt.subplot(121)
    plt.imshow(obj.image, interpolation='none')
    plt.contour(obj.image, colors='k')
    plt.subplot(122)
    plt.imshow(qim, interpolation='none')
    plt.contour(obj.image)
    #plt.contour(self.image_data[y-(Npx-1)/2:y+(Npx-1)/2+1,x-(Npx-1)/2:x+(Npx-1)/2+1])
    print x,y,x0,y0
    plt.savefig('tmp.png')
        '''


    x_img, y_img = obj.quad_image_pix_ref(x0, y0, x_images, y_images, fl)
    x_img -= (Npx-1)/2+1
    y_img -= (Npx-1)/2+1
    yg, xg =      np.mgrid[y-(Npx-1)/2:y+(Npx-1)/2+1,x-(Npx-1)/2:x+(Npx-1)/2+1]
    r, d = self.bigw.wcs_pix2world(xg-y0,yg-x0, 1)
    d = (d-dec)*3600.
    r = (r-ra) * np.cos(dec*np.pi/180.)*3600.
    fwhm = alpha*2.*np.sqrt(2.**(1/beta)-1.)*float(self.hdulist[0].header['PIXSCALE'])

    #vals = self.image_data[y-(Npx-1)/2:y+(Npx-1)/2+1,x-(Npx-1)/2:x+(Npx-1)/2+1]
    #vals = self.image_data[y-(Npx-1)/2:y+(Npx-1)/2+1,x-(Npx-1)/2:x+(Npx-1)/2+1]
    vals = obj.image
    min_ra = -4.
    max_ra = 4.
    min_dec = -4.
    max_dec = 4.
    vals -= nbkg
    vals*=ZP_flx
    #template = 0.*vals.copy()
    #for i in range(0,4):
    #    template[int(np.rint(y_im[i])-y+(Npx-1)/2),int(np.rint(x_im[i])-x+(Npx-1)/2)]=3.e-9
    vals0 = vals.copy()
    #print xg.shape, yg.shape, vals.shape
    ax1 = plt.subplot(3,2,2)
    #plt.pcolor(xg-x, yg-y, vals, cmap='jet', vmin=0., vmax=3.8e-9)
    plt.pcolor(r, d, vals, cmap='jet', vmin=0., vmax=3.8e-9)
    cb = plt.colorbar()
    plt.plot([min_ra,max_ra], [0.,0.], 'k--')
    plt.plot([0.,0.], [min_dec,max_dec], 'k--')

    plt.axis('equal')
    plt.xlim(max_ra, min_ra)
    plt.ylim(min_dec, max_dec)
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax1.yaxis.set_major_formatter(y_formatter)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax1.xaxis.set_major_formatter(x_formatter)
    #plt.xticks(rotation=30)
    #plt.title('x0 %1.3f y0 %1.3f'%(x0,y0))
    #plt.title('Data')
    plt.xlabel('$\Delta$RA,  arcseconds')
    plt.ylabel('$\Delta$Dec, arcseconds')
    ax1.text( 3.8,  3.5, 'Data', fontsize=30, color='white', fontweight='bold')
    ax1.text( 3.0,  0.0, 'M1', fontsize=30, color='white', fontweight='bold')
    ax1.text(-2.0, -0.5, 'M2', fontsize=30, color='white', fontweight='bold')
    ax1.text( 0.0,  2.0, 'S1', fontsize=30, color='white', fontweight='bold')
    ax1.text(0.5, -2.75, 'S2', fontsize=30, color='white', fontweight='bold')

    ax2 = plt.subplot(3,2,4)

    qim -= nbkg
    qim *= ZP_flx
    #plt.pcolor(xg-x,yg-y,qim, cmap='jet', vmin=0., vmax=3.8e-9)
    plt.pcolor(r,d,qim, cmap='jet', vmin=0., vmax=3.8e-9)
    cb = plt.colorbar()
    plt.plot([min_ra,max_ra], [0.,0.], 'k--')
    plt.plot([0.,0.], [min_dec,max_dec], 'k--')
    plt.axis('equal')
    plt.xlim(max_ra, min_ra)
    plt.ylim(min_dec, max_dec)
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax2.yaxis.set_major_formatter(y_formatter)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax2.xaxis.set_major_formatter(x_formatter)
    #plt.xticks(rotation=30)
    plt.xlabel('$\Delta$RA,  arcseconds')
    plt.ylabel('$\Delta$Dec, arcseconds')
    #plt.title('Model')
    ax2.text(3.8, 3.5, 'Model', fontsize=30, color='white', fontweight='bold')

    circle1=plt.Circle((3,-3),fwhm/2.,color='white', fill=False, linewidth=3, linestyle='dashed')
    plt.gcf().gca().add_artist(circle1)
    #plt.circles(3, -3, fwhm/2., alpha=1, lw=5, edgecolor='white', transform=ax2.transAxes)
    ax3 = plt.subplot(3,2,6)
    #plt.pcolor(xg-x,yg-y,vals-qim)
    plt.pcolor(r,d,vals-qim, vmin=-5.5e-10, vmax=5.5e-10)
    cb = plt.colorbar()
    plt.plot([min_ra,max_ra], [0.,0.], 'k--')
    plt.plot([0.,0.], [min_dec,max_dec], 'k--')
    plt.axis('equal')
    plt.xlim(max_ra, min_ra)
    plt.ylim(min_dec, max_dec)
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax3.yaxis.set_major_formatter(y_formatter)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax3.xaxis.set_major_formatter(x_formatter)
    #plt.xticks(rotation=30)
    plt.xlabel('$\Delta$RA,  arcseconds')
    plt.ylabel('$\Delta$Dec, arcseconds')
    ax3.text(3.8, 3.5, 'Residuals', fontsize=30, color='black', fontweight='bold')


################################################################
    
  def sub_stamp(self,ra, dec, pixels):
    x,y = self.bigw.wcs_world2pix(ra,dec,1)
    if(pixels%2==0):
	    #print 'even'
	    zoom_data = self.image_data[y-pixels/2:y+pixels/2,x-pixels/2:x+pixels/2]
    if(pixels%2==1):
	    #print 'odd'
	    zoom_data = self.image_data[y-(pixels-1)/2:y+(pixels-1)/2+1,x-(pixels-1)/2:x+(pixels-1)/2+1]
    return np.array(zoom_data)
    #plt.figure()

################################################################

  def get_exposure_time(self):
	  t1 = self.hdulist[0].header['UTSTART']
	  t2 = self.hdulist[0].header['UTSTOP']
	  #print t1+'000',t2
	  self.start_time = datetime.datetime.strptime(t1,"%H:%M:%S.%f")
	  self.end_time = datetime.datetime.strptime(t2,"%H:%M:%S.%f")
	  self.exposure_time = (self.end_time - self.start_time).seconds + (self.end_time - self.start_time).microseconds/1.e6
	  return self.exposure_time

# End FITManager Class Definition
################################################################
################################################################
################################################################

 
################################################################

def deVaucouleurs((x,y), _x0, _y0, r0, q = 1.0, posang=0.0):
    # Return a deVaucouleurs profile with integral normalized to unity.
    spa = np.sin(posang / 180.0 * np.pi)
    cpa = np.cos(posang / 180.0 * np.pi)

    # Define xprime coordinates in the rotated frame for convenience
    # This uses a standard rotation matrix
    xp = (x - _x0) * cpa - (y - _y0) * spa
    yp = (x - _x0) * spa + (y - _y0) * cpa
    # Define r with ellipticity and rotation
    r = np.sqrt(xp * xp / q / q + yp * yp)
    # The normalization was determined numerically. The normalized integral is good to a -1.4e-6 offset from unity.
    norm = 1./(q*22.665509867)
    return norm * np.exp(-7.669*( (r/r0)**0.25 - 1. ))

################################################################

def deVaucouleurs_model(_x0, _y0, F, r0, beta, shapelet_coeffs, nmax, mmax, q = 1.0, posang=0.0, npx=31, subsamp = 3, boundary_padding = 2):
    # This function produces a subpixelized and padded deVaucouleurs profile, convolves it with a Moffat PSF and performs the sub-pixel integration. 
    # The boundary is padded to reduce the boundary effects of the image convolution but the output has the correct image dimensions.
    # The steps are outlined below.
    #   1. Sub-pixelized and padded model of deVaucouleurs profile and moffat.
    #   2. Convolve the profiles
    #   3. Trim the boundaries.
    #   4. Sub-pixel integration
    #   5. Return output

    # 1. Sub-pixelized and padded model of deVaucouleurs profile and moffat.
    # these are the subpixel arrays we start with
    x = np.linspace(-npx / 2  + 0.5, npx / 2 + 0.5, npx * subsamp + 1  )
    y = np.linspace(-npx / 2  + 0.5, npx / 2 + 0.5, npx * subsamp + 1  )
    # this is how we pad them by a factor of two
    delta_x = np.diff(x)[0]
    delta_y = np.diff(y)[0]
    x = np.linspace( np.min(x) - delta_x*len(x)/2 , np.max(x) + delta_x*len(x)/2, len(x)*2  )
    y = np.linspace( np.min(y) - delta_y*len(y)/2 , np.max(y) + delta_y*len(y)/2, len(y)*2  )

    # define the grid coordinates for the subsampled model
    xg, yg = np.meshgrid(x, y)

    # compute the model
    img = deVaucouleurs((xg,yg), _x0, _y0, r0, q=q, posang=posang) 

    # get the shapelet model and normalize it to unity, including the subsampling
    PSF_model = psf.shapelet_image(xg, yg, nmax, mmax, 0., 0., beta,  shapelet_coeffs)
    PSF_model /= subsamp**2 * 2.*beta*np.sqrt(np.pi)

    # 2. Convolve the profiles
    conv = np.fft.irfft2(np.fft.rfft2(PSF_model) * np.fft.rfft2(img)) # FFT2 is much faster than scipy convolve and gives the same results.
    conv = np.fft.fftshift(conv)

    # 3. Trim the boundaries.
    trm = conv[ len(x)/4: -len(x)/4, len(y)/4: -len(y)/4 ]

    # 4. Sub-pixel integration (lifted from Curtis' code. See psf.moffat_kernel function for more detailed comments.)    
    x2 = np.linspace(-npx / 2  + 0.5, npx / 2 + 0.5, npx * subsamp + 1  )
    y2 = np.linspace(-npx / 2  + 0.5, npx / 2 + 0.5, npx * subsamp + 1  )
    x2d, y2d = np.meshgrid(x2, y2)
    subkern=trm.copy()
    result4d = np.zeros((npx, npx, subsamp + 1, subsamp + 1))
    kern4d = integrate.make4d(subkern, npx, npx, subsamp)
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    subsamp_img = integrate.simp4d(kern4d, dx, dy)

    # rotate and flip to be consistent with the star and quasar image models
    subsamp_img = np.rot90(subsamp_img)
    subsamp_img = np.flipud(subsamp_img)
    subsamp_img *= F # This value of the flux should be in units of integrated sum_CCD

    '''
    print 'beta',beta
    print 'np.sum(PSF_model)',np.sum(PSF_model)
    print 'np.sum(img)',np.sum(img)
    print 'np.sum(conv)',np.sum(conv)
    print 'np.sum(trm)',np.sum(trm)
    print 'np.sum(integrate.simp4d(kern4d, dx, dy))', np.sum(integrate.simp4d(kern4d, dx, dy))
    plt.figure()
    plt.subplot(221)
    plt.imshow(img, interpolation='none', cmap='viridis')
    plt.colorbar()
    plt.subplot(222)
    plt.imshow(conv, interpolation='none', cmap='viridis')
    plt.colorbar()
    plt.subplot(223)
    plt.imshow(trm, interpolation='none', cmap='viridis')
    plt.colorbar()
    plt.subplot(224)
    plt.imshow(subsamp_img, interpolation='none', cmap='viridis')
    plt.colorbar()
    plt.savefig('psf_model.png')
    #exit()
    '''
    return subsamp_img
      
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

class SourceImage:
  def __init__(self, FM, ra, dec, nmax, mmax, pixels):
    #self.FM = FM
    # you will set the relevant parameters as class members here, not all of FM.
    self.image  = FM.sub_stamp(ra, dec, pixels)
    self.gain   = FM.hdulist[0].header['GAIN']
    self.fits_file_name = FM.fits_file_name
    self.ra     = ra
    self.dec    = dec
    self.nmax   = nmax
    self.mmax   = mmax
    self.pixels = pixels


    # Mask bad values
    # curve fit does not respond to masked arrays. The data has to be fed with these points actually excluded.
    # this is a somewhat time consuming change to make and validate.
    # self.image = np.ma.masked_where(self.image<=0,self.image)
    '''
    # find ratio of nearest neighbors
    Nx,Ny = self.image.shape
    mask = np.zeros(self.image.shape)
    for i in range(0,Nx):
      for j in range(0,Ny):
        min_ratio=1.e9
        for di in [-1,0,1]:
          for dj in [-1,0,1]:
            if(i+di>=0 and i+di<Nx and j+dj>0 and j+dj<Nx and (di!=0 and dj!=0)):
		#print i,j,di,dj 		
		ratio = self.image[i][j]/self.image[i+di][j+dj]
	        #print 'ratio', ratio
		if(ratio<min_ratio):
			min_ratio = ratio
	#if(min_ratio>1.): print min_ratio
        if(min_ratio>1.5):
		print 'MASKING PIXEL',i,j
		mask[i][j]=1
    self.image = np.ma.masked_where(mask==1,self.image)
    '''
  '''
  def stacked_twoD_Moffat_chi(self, theta, readnoise):
      chisq=0.
      for k in range(0,len(self.PSFimage_data))
  '''

  #############################################################################

  def fit_shapelets(self, readnoise, nmax, mmax, verbose=False, display=0):
    #if(verbose): print 'fitting shapelets'

    x=np.arange(0.,len(self.image))
    y=np.arange(0.,len(self.image))
    self.xg, self.yg = np.mgrid[:len(self.image), :len(self.image)]

    # estimate centroid
    self.x_max, self.y_max = np.unravel_index(self.image.argmax(), self.image.shape)
    #print '\tx_max, y_max', self.x_max, self.y_max
    x0_guess = self.x_max
    y0_guess = self.y_max
    x0_guess = float(x0_guess - (len(x)- 1) / 2) 
    y0_guess = float(y0_guess - (len(y)- 1) / 2)
    #print '\tx0_guess, y0_guess', x0_guess, y0_guess

    # estimate background counts
    background_count_guess = np.min(self.image)
    #print '\tbackground_count_guess ', background_count_guess  

    # estimate fwhm
    pk = np.max(self.image)
    fwhm=1.
    #print '\tpk, fwhm', pk, fwhm
    for dx in range(0,15):
        val_a = self.image[self.x_max + dx, self.y_max]
        val_b = self.image[self.x_max + dx+1, self.y_max]
        if(val_a>=pk/2. and val_b<=pk/2.):
	     fwhm = dx+0.5
	     break
    #print '\tcrude fwhm', fwhm

    #if(verbose): print 'initial guess', self.initial_guess_moffat
    #if(verbose): print len(self.image)

    # guess shapelet beta parameter
    guess_beta = fwhm/1.4
    # initizalize shapelet parameter array
    #print 'NMAX', nmax
    #params = np.zeros(((nmax+1)*(nmax+2))//2) 
    params = np.zeros(psf.num_params(nmax, mmax))
    params[0]=self.image[self.x_max,self.y_max]*np.sqrt(np.pi)*guess_beta # peak value Gaussian
    #print 'NMAX', nmax, len(params)
    #print params
    #xo = float(x0_guess - (len(x)- 1) / 2)
    #yo = float(y0_guess - (len(y)- 1) / 2)
    self.initial_guess_shapelets = np.concatenate([[x0_guess, y0_guess, guess_beta, background_count_guess],params])

    #print self.initial_guess_shapelets
    #print 'xg.shape', self.xg.shape
    m = psf.shapelet_image(self.xg - (len(x)- 1) / 2, self.yg - (len(y)- 1) / 2, nmax, mmax, x0_guess, y0_guess, guess_beta, params) + background_count_guess 
    d = m-self.image
    spim = self.polarShapelet((self.xg, self.yg), *self.initial_guess_shapelets).reshape(self.xg.shape)
    dspim = spim-self.image
    #gain = self.FM.hdulist[0].header['GAIN']
    sigmas = np.sqrt( readnoise**2 + self.image/self.gain )

    self.fit_ok = True
    #popt, pcov = opt.curve_fit(self.polarShapelet, (self.xg, self.yg), self.image.ravel(), p0=self.initial_guess_shapelets)
    popt, pcov = opt.curve_fit(self.polarShapelet, (self.xg, self.yg), self.image.ravel(), sigma = sigmas.ravel(), absolute_sigma=True, p0=self.initial_guess_shapelets)
    #popt, pcov = opt.curve_fit(self.polarShapelet, (self.xg, self.yg), self.image.ravel(), sigma = sigmas.ravel(), p0=popt)
    if(np.isnan(pcov).any() or np.isinf(pcov).any()):
        self.fit_ok = False
    #print popt
    self.shapelet_fit_parameters = popt
    self.shapelet_fit_covariance = pcov
    #spim_ss = psf.shapelet_kernel(popt, nmax, mmax, len(x), len(y))
    spim_ss = self.polarShapelet((self.xg, self.yg), *popt).reshape(self.xg.shape)
    dspim_ss = spim_ss-self.image
    self.shapelet_fit_image = spim_ss
    chi = dspim_ss/sigmas
    chi_sq =  np.sum( chi**2 )/float(len(sigmas.ravel()) - len(params))
    max_chi = np.max( np.abs( chi ) )
    #print 'fitted: chis_sq, max_chi', chi_sq, max_chi  

    #print 'chi from mad', 1.4826*np.median(np.abs(chi - np.median(chi)))
    #print '\tpopt',popt
    #print '\tpcov',np.sqrt(pcov[0][0]), np.sqrt(pcov[0][0])/popt[0], np.linalg.det(pcov)
    #EXTRACT CCD COUNT INTEGRAL FROM PSF COEFFICIENTS
    # get a list of coefficients with m=0
    nm_list = psf.get_nm(nmax, mmax)
    fn0_list = []
    fn_indices = []
    beta_fit = popt[2]
    for k in range(0,len(nm_list)):
        if(nm_list[k][1]==0):
            fn0_list.append(popt[k+4])
            fn_indices.append(k+4)
    #print 'len(fn0_list), nmax',len(fn0_list), nmax     
    S_CCD = np.sum(fn0_list)*2*np.sqrt(np.pi)*beta_fit
    #print 'pcov',pcov
    #print 'pcov[2][fn_indices]', pcov[2][fn_indices]
    #print 'fn_indices', fn_indices
    #print 'pcov.shape', pcov.shape
    ixgrid = np.ix_(fn_indices, fn_indices)
    #print 'ixgrid', ixgrid
    #print 'pcov[ixgrid]',pcov[ixgrid], pcov[ixgrid].shape

    #print '(1)', 4. * np.pi * np.sum(fn0_list)**2 * pcov[2][2]  
    #print '(2)', 4. * np.pi * popt[2] * np.sum( np.outer( pcov[2][fn_indices] , popt[fn_indices] ) )
    #print '(3)', 4. * np.pi * popt[2]**2 * np.sum( pcov[ np.ix_( fn_indices, fn_indices ) ] )

    sig2_S_CCD  = 4. * np.pi * np.sum(fn0_list)**2 * pcov[2][2]
    sig2_S_CCD += 4. * np.pi * popt[2] * np.sum( np.outer( pcov[2][fn_indices] , popt[fn_indices] ) )
    sig2_S_CCD += 4. * np.pi * popt[2]**2 * np.sum( pcov[ np.ix_( fn_indices, fn_indices ) ] )

    #print 'S_CCD', S_CCD, np.sum(spim_ss - popt[3])*1.266
    #print 'sig2_S_CCD', sig2_S_CCD, np.sqrt(sig2_S_CCD), np.sqrt(sig2_S_CCD)/S_CCD 
    # add some conditions to identify a bad fit here

    if(display>2):
        plt.figure(figsize=(12,6))
        plt.subplot(121)
        plt.imshow(m, interpolation='none', cmap='viridis')
        plt.colorbar()
        plt.title('Guess Image')
        plt.subplot(122)
        vv = np.max(np.abs(d))
        plt.imshow(d, interpolation='none', vmin = -vv, vmax=vv, cmap='seismic')
        plt.colorbar()
        plt.xlabel('pixel')
        plt.ylabel('pixel')
        plt.title('Guess - Data')

        plt.figure(figsize=(12,6))
        plt.subplot(121)
        plt.imshow(spim, interpolation='none', cmap='viridis')
        plt.colorbar()
        plt.title('Guess Image')
        plt.subplot(122)
        vv = np.max(np.abs(dspim))
        plt.imshow(dspim, interpolation='none', vmin = -vv, vmax=vv, cmap='seismic')
        plt.colorbar()
        plt.xlabel('pixel')
        plt.ylabel('pixel')
        plt.title('Guess - Data')

        plt.figure(figsize=(12,6))
        plt.subplot(121)
        plt.imshow(spim_ss, interpolation='none', cmap='viridis')
        plt.colorbar()
        plt.title('Fitted Image')
        plt.subplot(122)
        vv = np.max(np.abs(dspim_ss))
        plt.imshow(dspim_ss, interpolation='none', vmin = -vv, vmax=vv, cmap='seismic')
        plt.colorbar()
        plt.xlabel('pixel')
        plt.ylabel('pixel')
        plt.title('Fit - Data')

        plt.figure(figsize=(12,6))
        plt.subplot(121)
        plt.hist((dspim_ss/sigmas).ravel())
        plt.title('Chi Values')
        plt.subplot(122)
        vv = np.max(np.abs(dspim_ss/sigmas))
        plt.imshow(dspim_ss/sigmas, interpolation='none', vmin = -vv, vmax=vv, cmap='seismic')
        plt.colorbar()
        plt.xlabel('pixel')
        plt.ylabel('pixel')
        plt.title('Chi Values')

        plt.figure(figsize=(12,6))
        psd2 = np.abs(np.fft.fft2(dspim_ss/sigmas))  
        plt.subplot(121)
        plt.hist(psd2.ravel())
        plt.title('Residuals Power Spectrum ')
        plt.subplot(122)
        vv = np.max(np.abs(psd2))
        plt.imshow(np.fft.fftshift(psd2), interpolation='none', cmap='viridis', norm=LogNorm())
        plt.colorbar()
        plt.xlabel('pixel')
        plt.ylabel('pixel')
        plt.title('Power Spectrum of Residuals')

    return S_CCD, np.sqrt(sig2_S_CCD), chi_sq, max_chi, popt

  #############################################################################

  def polarShapelet(self, (x, y), *parameters):
    m = psf.shapelet_kernel(parameters, self.nmax, self.mmax, nx = len(x), ny = len(y))
    bkg = parameters[3]
    m = np.rot90(m)
    m = np.flipud(m)
    #m = amplitude*PSF_model + offset
    if(parameters[2]<0. or parameters[3]<0.): 
        m+=1.e9
    return m.ravel() + bkg


  #############################################################################


  def quad_image_ref(self, ra_offset, dec_offset):
	# values from Kochanek 2006
	ra  = np.array([0.,  2.467, 1.476,  0.939])
	dec = np.array([0., -0.603, 0.553, -1.614])
	ra*=-1.
	return ra-ra_offset, dec-dec_offset

  #############################################################################

  def quad_image_pix_ref(self, x0, y0, x_image, y_image, flip):
    x1 = y_image[0]+x0
    x2 = y_image[1]+x0
    x3 = y_image[2]+x0
    x4 = y_image[3]+x0
    y1 = x_image[0]+y0
    y2 = x_image[1]+y0
    y3 = x_image[2]+y0
    y4 = x_image[3]+y0

    x_vals = np.array([x1,x2,x3,x4]) 
    y_vals = np.array([y1,y2,y3,y4])
    #if(flip):
    #	y_vals = 15.-1.*(y_vals-15.)
    return x_vals, y_vals

  #############################################################################

  def quad_image_model(self, (xg,yg, beta, shapelet_coeffs, N_pix, flip), x0, y0, _x1, _x2, _x3, _x4, _y1, _y2, _y3, _y4, amp0, amp1, amp2, amp3, N_bkg, _xlens=0., _ylens=0., amp_lens=0., r0_lens=1., q_lens=1., posang_lens=0.):
    #t0 = time.clock()
    # image x's and y's are flipped. unflip them and add the centroid offset
    x1 = _y1 + x0
    x2 = _y2 + x0
    x3 = _y3 + x0
    x4 = _y4 + x0
    y1 = _x1 + y0
    y2 = _x2 + y0
    y3 = _x3 + y0
    y4 = _x4 + y0
    xlens = _ylens + x0
    ylens = _xlens + y0
    '''
    plt.figure()
    plt.subplot(221)
    plt.plot([_x1],[_y1], 'o')
    plt.plot([_x2],[_y2], 'o')
    plt.plot([_x3],[_y3], 'o')
    plt.plot([_x4],[_y4], 'o')
    plt.plot([_xlens],[_ylens], 'o')
    plt.subplot(222)
    plt.plot([x1],[y1], 'o')
    plt.plot([x2],[y2], 'o')
    plt.plot([x3],[y3], 'o')
    plt.plot([x4],[y4], 'o')
    plt.plot([xlens],[ylens], 'o')
    plt.savefig('positions_in_qim_%1.5f.png'%x0)
    '''
   
    #print 'x0,y0', x0,y0
    #print 'x1,y1', x1,y1
    #print 'x2,y2', x2,y2
    #print 'x3,y3', x3,y3
    #print 'x4,y4', x4,y4
   
    # NOTE: THIS FUNCTION IS DEFINED TWICE. INTERNALLY IN EACH CASE BUT WE SHOULD IMPROVE THIS BY FINDING A WAY TO DEFINE IT ONCE
    def shapeletMagnitude((x, y), _x0, _y0, bkg, sum_CCD):
      #print 'sum_CCD', sum_CCD
      parameters = np.concatenate([[_x0,_y0, beta, bkg], sum_CCD*shapelet_coeffs])
      m = psf.shapelet_kernel(parameters, self.nmax, self.mmax, nx = len(x), ny = len(y))
      bkg = parameters[3]
      m = np.rot90(m)
      m = np.flipud(m)
      #m = amplitude*PSF_model + offset
      if(parameters[2]<0. or parameters[3]<0.): 
          m+=1.e9
      return m.ravel() + bkg

    #p0  =  twoD_Moffat((xg, yg), amp0, alpha, beta, x1, y1, 0)
    #p1  =  twoD_Moffat((xg, yg), amp1, alpha, beta, x2, y2, 0)
    #p2  =  twoD_Moffat((xg, yg), amp2, alpha, beta, x3, y3, 0)
    #p3  =  twoD_Moffat((xg, yg), amp3, alpha, beta, x4, y4, 0)

    p0  =  shapeletMagnitude((xg, yg), x1, y1, 0., amp0)
    p1  =  shapeletMagnitude((xg, yg), x2, y2, 0., amp1)
    p2  =  shapeletMagnitude((xg, yg), x3, y3, 0., amp2)
    p3  =  shapeletMagnitude((xg, yg), x4, y4, 0., amp3)
    #print 'max p0', np.max(p0)
    '''
    plt.figure()
    plt.subplot(221)
    plt.imshow(p0.reshape(N_pix, N_pix), interpolation='none')
    plt.colorbar()
    plt.subplot(222)
    plt.imshow(p1.reshape(N_pix, N_pix), interpolation='none')
    plt.colorbar()
    plt.subplot(223)
    plt.imshow(p2.reshape(N_pix, N_pix), interpolation='none')
    plt.colorbar()
    plt.subplot(224)
    plt.imshow(g, interpolation='none')
    plt.colorbar()
    plt.savefig('test_images.png')
    plt.show()
    exit()
    '''
    if(amp_lens!=0.):
        # print 'x_lens, y_lens, amp_lens, r0_lens, q_lens, posang_lens', x_lens, y_lens, amp_lens, r0_lens, q_lens, posang_lens
        # g   = deVaucouleurs_model(xlens, ylens, amp_lens, r0_lens, alpha, beta, q = q_lens, posang=posang_lens, npx=31)
        # print 'in quad_imag_model, xlens, ylens', xlens, ylens
        g   = deVaucouleurs_model(xlens, ylens, amp_lens, r0_lens, beta, shapelet_coeffs, self.nmax, self.mmax, q = q_lens, posang=posang_lens, npx=31)
        #print 'g.shape', g.shape, 'p0.shape', p0.shape 
        '''
        plt.figure()
        x1,y1 = np.flipud(np.rot90([[_x1,_y1]]))
        x2,y2 = np.flipud(np.rot90([[_x2,_y2]]))
        x3,y3 = np.flipud(np.rot90([[_x3,_y3]]))
        x4,y4 = np.flipud(np.rot90([[_x4,_y4]]))
        xlens,ylens = np.flipud(np.rot90([[_xlens,_ylens]]))

        plt.plot([x1+(len(xg)-1)/2],[y1+(len(yg)-1)/2], 'o')
        plt.plot([x2+(len(xg)-1)/2],[y2+(len(yg)-1)/2], 'o')
        plt.plot([x3+(len(xg)-1)/2],[y3+(len(yg)-1)/2], 'o')
        plt.plot([x4+(len(xg)-1)/2],[y4+(len(yg)-1)/2], 'o')
        plt.plot([xlens+(len(xg)-1)/2],[ylens+(len(yg)-1)/2], 'o')
        if np.max(p0) > 0. and np.max(p1) > 0. and np.max(p2) > 0. and np.max(p3) > 0. :
          plt.contour(p0.reshape(31,31), levels=np.linspace(0.5*np.max(p0), np.max(p0), 10))
          plt.contour(p1.reshape(31,31), levels=np.linspace(0.5*np.max(p1), np.max(p1), 10))
          plt.contour(p2.reshape(31,31), levels=np.linspace(0.5*np.max(p2), np.max(p2), 10))
          plt.contour(p3.reshape(31,31), levels=np.linspace(0.5*np.max(p3), np.max(p3), 10))
          plt.contour(g, levels=np.linspace(0.5*np.max(g), np.max(g), 10))
          plt.colorbar()
          plt.savefig('gal_in_qim.png')
        '''
        return np.ravel(g) + (p0+p1+p2+p3) + N_bkg
    # otherwise return quasar images only
    return (p0+p1+p2+p3) + N_bkg

# End SourceImage Class Definition
################################################################
################################################################
################################################################


################################################################

def gauss_1d(x, amp, sig, x0):
	return amp*np.exp(-(x-x0)**2/2./sig**2)

################################################################

def magnitude(value):
  return -2.5*np.log10(value)

################################################################

####################################################################
## PHOTOMETRY FUNCTIONS ############################################
####################################################################

def estimate_total_light(obj, N_bkg, sigma_read, display=0, out='out'):
  estimate, uncertainty, chi_sq, max_chi, fit_parms = obj.fit_shapelets(sigma_read, obj.nmax, obj.mmax, display=display)
  res = obj.image-obj.shapelet_fit_image
  chi = res/(np.sqrt(obj.image + sigma_read**2))
  if(display >= 3):
	fig = plt.figure(figsize=(8,11))
	plt.subplot(321)
	plt.imshow(obj.image, interpolation='none', vmin=np.min(obj.image), vmax=np.max(obj.image), cmap='viridis')
	plt.colorbar()
	plt.title('data')
	plt.subplot(322)
	plt.imshow(obj.shapelet_fit_image, interpolation='none', vmin=np.min(obj.image), vmax=np.max(obj.image), cmap='viridis')
	plt.colorbar()
	plt.title('Shapelet Fit\n(n_max=%d, m_max=%d)'%(obj.nmax, obj.mmax))
	plt.subplot(323)
	plt.imshow(res, interpolation='none', vmin=-np.max(np.max(res)),  vmax=np.max(np.max(res)), cmap='seismic')
	plt.colorbar()
	plt.title('residuals')
	plt.subplot(324)
	plt.imshow(chi, interpolation='none', vmin=-np.max(np.abs(chi)), vmax=np.max(np.abs(chi)), cmap='seismic')
	mx = np.max(np.abs(chi))
	mn = -mx
	plt.colorbar()
	plt.title('chi values')
	plt.subplot(313)
	plt.hist(np.ravel(chi), bins=30)
	plt.xlim(mn,mx)
	plt.suptitle(obj.fits_file_name.split('/')[-1]+'\n'+out, fontsize=20)
	plt.xlabel('chi values')
	plt.subplots_adjust(bottom=0.05)
	plt.title('chi distribution')
	plt.savefig(out+'.png', dpi=50)
	plt.close(fig)
  return estimate, uncertainty, chi_sq, max_chi

################################################################
def estimate_PSF(FM, PSF_star_table_filename, sigma_read, nmax, mmax, npx, display=0, out = 'out'):
  # PSF stars
  PSF_star_table = ascii.read(PSF_star_table_filename)

  # KEEP TRACK OF WHICH BAND THE IMAGES ARE AND USE THE CORRESPONDING APASS REFERENCE VALUES.
  print '\nestimate_PSF'
  filt = FM.hdulist[0].header['FILTER'] # gp or rp
  print '\tPSF list file :', PSF_star_table_filename
  print '\tFilter :', filt
  num_psf_stars = len(PSF_star_table)
  num_psf_parms = 4+psf.num_params(nmax, mmax) # additional 4 are x0,y0, beta, background counts
  # INITIALIZE ARRAYS FOR THE PSF STAR MAGNITUDES AND PSF PARAMETERS
  beta             = np.zeros(num_psf_stars)
  beta_error       = np.zeros(num_psf_stars)
  bkg              = np.zeros(num_psf_stars)
  bkg_error        = np.zeros(num_psf_stars)
  fit_parms        = np.zeros((num_psf_stars, num_psf_parms))
  fit_cov          = np.zeros((num_psf_stars, num_psf_parms, num_psf_parms))
  S_CCD_list       = np.zeros(num_psf_stars)
  sig_S_CCD_list   = np.zeros(num_psf_stars)
  chi_sq_list      = np.zeros(num_psf_stars)
  max_chi_list     = np.zeros(num_psf_stars)
  fit_success      = np.zeros(num_psf_stars, dtype=np.bool)

  for k in range(0,len(PSF_star_table)):
	print '\tPSF star index = %d'%k
	print '\t\tra , dec'.ljust(15), '= %1.5f , %1.5f'%(PSF_star_table['radeg'][k], PSF_star_table['decdeg'][k])

	obj = SourceImage(FM, PSF_star_table['radeg'][k], PSF_star_table['decdeg'][k], nmax, mmax, npx)
	N_bkg = np.median(obj.image)
	outTag = out+'_PSF_src_%d'%k
	try:
		intg, intg_unc, chi_sq, max_chi = estimate_total_light(obj, N_bkg, sigma_read, display=display, out=outTag)
	except:
		print '\t\tFIT FAILED' 

	if(obj.fit_ok):
	  print '\t\tFIT SUCCEEDED'
	  print '\t\tFlux'.ljust(15), '= %1.0f +/- %1.0f'%(intg, intg_unc)
	  print '\t\tchi^2'.ljust(15), '= %1.2f'%chi_sq
	  print '\t\tmax|chi|'.ljust(15), '= %1.2f'%max_chi
	  S_CCD_list[k]     = intg
	  sig_S_CCD_list[k] = intg_unc
	  #beta[k]          = obj.shapelet_fit_parameters[2]
	  #psf_parms[k]     = obj.shapelet_fit_parameters[4:]
	  #psf_parm_err[k]  = np.sqrt(np.diagonal(obj.shapelet_fit_covariance)[4:])
	  fit_parms[k]      = obj.shapelet_fit_parameters
	  fit_cov[k]        = obj.shapelet_fit_covariance
	  chi_sq_list[k]    = chi_sq
	  max_chi_list[k]   = max_chi

	  fit_success[k] = True

  print '\tfit_success', fit_success

  #beta_med  = np.median(beta)
  #beta_unc  = 1.4826*np.median(np.abs(beta -  beta_med))

  # SAVE OUTPUTS TO NPZ FILE
  # print '\tBeta'.ljust(15), '%1.3f +/- %1.3f'%(beta_med,  beta_unc)
  npz_out = out + '_PSF.npz'
  np.savez( npz_out,
            imageFile       = FM.fits_file_name,
            outFileTag      = out,
            mjd_obs         = float(FM.hdulist[0].header['MJD-OBS']),
            readnoise       = FM.readnoise,
            PSF_list_file   = PSF_star_table_filename,
            fit_parms       = fit_parms,
            fit_cov         = fit_cov,
            nmax            = nmax,
            mmax            = mmax,
            filter          = filt,
            S_CCD           = S_CCD_list,
            sig_S_CCD       = sig_S_CCD_list,
            chi_sq          = chi_sq_list, 
            max_chi         = max_chi_list,
            fit_success     = fit_success
  )

  return 0

def Get_Median_PSF_Parameters(PSF_file, chi_sq_cut, max_chi_cut, nmax, mmax ):
    # load file
    f = np.load(PSF_file)
    # set cuts
    cut = np.logical_and(f['chi_sq'] < chi_sq_cut, f['max_chi'] < max_chi_cut) 
    cut = np.logical_and(f['fit_success'], cut)

    # fill normalized parameters
    normed_parms = []
    normed_parm_err = []
    good_count = 0
    nm_list = np.array(psf.get_nm(nmax, mmax))
    #print 'nm_list', nm_list
    #print 'nm_list[:,1]', nm_list[:,1]
    #print 'nm_list.shape', nm_list.shape

    fn0_indices = np.where(nm_list[:,1]==0)[0]
    #print 'fn0_indices', fn0_indices
    for k in range(0,len(cut)):
        if(cut[k]):
            good_count += 1
            normed_parms.append( f['fit_parms'][k,4:] / f['S_CCD'][k] ) # * 2.*np.sqrt(np.pi)* f['fit_parms'][k,2])
            normed_parm_err.append( np.sqrt(np.diagonal(f['fit_cov'][k])[4:]) / f['S_CCD'][k])# * 2.*np.sqrt(np.pi)* f['fit_parms'][k,2] )
	    #print ''
            #print 'k, S_CCD', k, f['S_CCD'][k]
	    #print len(normed_parms[-1])
	    #print 'f[\'fit_parms\'][k,:]', f['fit_parms'][k,:]
	    #print 'fn0_indices', fn0_indices
            #print f['fit_parms'][k,4+fn0_indices]/f['S_CCD'][k] - normed_parms[-1][fn0_indices]
            #print 'np.sum(f[\'fit_parms\'][k,4+fn0_indices] / f[\'S_CCD\'][k])*2.*np.sqrt(np.pi)* f[\'fit_parms\'][k,2]',np.sum(f['fit_parms'][k,4+fn0_indices]/f['S_CCD'][k])*2.*np.sqrt(np.pi)* f['fit_parms'][k,2], f['S_CCD'][k], np.sum(f['fit_parms'][k,4+fn0_indices])*2.*np.sqrt(np.pi)* f['fit_parms'][k,2]/f['S_CCD'][k]
	    #print 'np.sum(f[\'fit_parms\'][k,4+fn0_indices])*2.*np.sqrt(np.pi)* f[\'fit_parms\'][k,2]',np.sum(f['fit_parms'][k,4+fn0_indices])*2.*np.sqrt(np.pi)* f['fit_parms'][k,2], f['S_CCD'][k], np.sum(f['fit_parms'][k,4+fn0_indices])*2.*np.sqrt(np.pi)* f['fit_parms'][k,2]/f['S_CCD'][k]
            #print 'normed_parms[-1]', normed_parms[-1]
	    #print 'np.sum(normed_parms[-1][fn0_indices])*2.*np.sqrt(np.pi)* f[\'fit_parms\'][k,2]',np.sum(normed_parms[-1][fn0_indices])*2.*np.sqrt(np.pi)* f['fit_parms'][k,2]

    print 'Get_Median_PSF_Parameters'
    print '\tNumber of Good PSF Fits', good_count
    normed_parms = np.array(normed_parms)
    normed_parm_err = np.array(normed_parm_err)
    PSF_parm_med = np.median(normed_parms, axis=0)
    PSF_parm_MAD = 1.4826* np.median ( np.abs ( normed_parms - np.median(normed_parms, axis=0)), axis=0)

    # we want to enforce sum(f_n0)*2*sqrt(pi)*beta = 1, so we don't take medians or means for the final beta.
    #beta_med = np.median(f['fit_parms'][:,2][cut] )
    #beta_MAD = 1.4826*np.median(np.abs ( f['fit_parms'][:,2][cut] - np.median(f['fit_parms'][:,2][cut])) )
    f.close()

    # rather we set
    beta = 1./( np.sum(PSF_parm_med[fn0_indices])*2.*np.sqrt(np.pi) )
    beta_unc = np.sqrt(np.sum(PSF_parm_MAD[fn0_indices]**2)) / ( np.sum(PSF_parm_med[fn0_indices])**2 *2.*np.sqrt(np.pi) ) # not totally  proper error propagation but it does give and estimate
    #print 'beta_med = %1.2f +/- %1.2f'%(beta_med, beta_MAD) 
    print '\tbeta     = %1.2f +/- %1.2f'%(beta, beta_unc)
    #print 'np.sum(PSF_parm_med[fn0_indices])*2.*np.sqrt(np.pi)*beta_med', np.sum(PSF_parm_med[fn0_indices])*2.*np.sqrt(np.pi)*beta
    #print PSF_parm_med
    #print PSF_parm_MAD
    return beta, beta_unc, PSF_parm_med, PSF_parm_MAD, good_count


def APASS_zero_points(FM, APASS_table, APASS_rejects, sigma_read, display=0, out = 'out'):
  # KEEP TRACK OF WHICH BAND THE IMAGES ARE AND USE THE CORRESPONDING APASS REFERENCE VALUES.
  print '\nAPASS_zero_points'
  filt = 'Sloan_r'
  filt_err = 'r_err'
  if( FM.hdulist[0].header['FILTER'] == 'gp'):
    filt = 'Sloan_g'
    filt_err = 'gerr'
  # MAKE LISTST FOR THE CALIBRATED APASS MAGNITUDES AND THE INSTRUMENT MAGNITUDES
  m_APASS = []
  m_I     = []
  APASS_index_list = []
  m_APASS_unc = []
  m_I_unc     = []
  beta  = []
  beta_error  = []
  bkg=[]
  bkg_error=[]
  psf_parms = []
  psf_parm_err = []
  S_CCD_list = []
  sig_S_CCD_list = []
  chi_sq_list = []
  max_chi_list = []
  pass_quality_cut = []
  for k in range(0,len(APASS_table)):
	#print ''
	#print APASS_table[filt][k], APASS_table[filt_err][k], k 
	#print APASS_table[filt][k]=='NA', APASS_table[filt_err][k]=='0', k in APASS_rejects
	#print filt
	if(APASS_table[filt][k]=='NA' or APASS_table[filt_err][k]=='0' or k in APASS_rejects):
		continue
	print '\tAPASS index = %d'%k
	print '\t\tra , dec'.ljust(15), '= %1.5f , %1.5f'%(APASS_table['radeg'][k], APASS_table['decdeg'][k])
	print ('\t\t%s mag'%(filt)).ljust(15), '= %s'%(APASS_table[filt][k])
	obj = SourceImage(FM, APASS_table['radeg'][k], APASS_table['decdeg'][k], 31)
	N_bkg = np.median(obj.image)
	outTag = out+'_src_%d'%k
	#plt.show()
	#'''
	try:
		intg, intg_unc, chi_sq, max_chi = estimate_total_light(obj, N_bkg, sigma_read, display=display, out=outTag)
	except:
		print '\t\tFIT FAILED' 
		continue

	if(obj.fit_ok):
	  print '\t\tFIT SUCCEEDED'
	  print '\t\tFlux'.ljust(15), '= %1.0f +/- %1.0f'%(intg, intg_unc)
	  print '\t\tchi^2'.ljust(15), '= %1.2f'%chi_sq
	  print '\t\tmax|chi|'.ljust(15), '= %1.2f'%max_chi
	  APASS_index_list.append(k)
	  S_CCD_list.append(intg)
	  sig_S_CCD_list.append(intg_unc)
	  #m_I.append(flux2magnitude(intg))
	  #m_I_unc.append(fluxErr2magErr(intg, intg_unc))
	  #m_APASS.append(float(APASS_table[filt][k]))
	  #m_APASS_unc.append(float(APASS_table[filt_err][k]))
	  beta.append(obj.shapelet_fit_parameters[2])
	  #beta_error.append(np.sqrt(obj.shapelet_fit_covariance[2][2]))
	  #bkg.append(obj.shapelet_fit_parameters[3])
	  #bkg_error.append(np.sqrt(obj.shapelet_fit_covariance[3][3]))
	  psf_parms.append(obj.shapelet_fit_parameters[4:])
	  psf_parm_err.append(np.sqrt(np.diagonal(obj.shapelet_fit_covariance)[4:]))
	  chi_sq_list.append(chi_sq)
	  max_chi_list.append(max_chi)

	  pass_quality_cut.append(True)
	  #'''
	  # DO NOT USE FOR PSF PARAMETER ESTIMATION UNLESS IT PASSES THESE QUALITY CUTS
	  chi_sq_cut  = 1.3  * 1.
	  max_chi_cut = 4.5  * 1.
	  if(chi_sq > chi_sq_cut or max_chi > max_chi_cut):
	    print '\t\tQUALITY CUTS FAILED'
	    print '\t\t\tchi^2    %1.2f > %1.2f'%(chi_sq,  chi_sq_cut)
	    print '\t\t\tmax|chi| %1.2f > %1.2f'%(max_chi, max_chi_cut)
	    pass_quality_cut[-1] = False
	    #continue

  '''
  #plt.show()
  # ESTIMATE ZERO POINTS FOR EACH APASS SOURCE
  ZP = np.array(m_APASS)-np.array(m_I)
  # USE ESTIMATED UNCERTINTIES TO WEIGHT THE MEAN AND RMS OF ZERO POINTS
  weight = 1./np.sqrt(np.array(m_I_unc)**2 + np.array(m_APASS_unc)**2)
  # ESTIMATE THE MEAN AND WEIGHTED RMS OF THE ZERO POINTS 
  #ZP_mean = sum(ZP)/len(ZP)
  ZP_rms = np.sqrt(sum(ZP**2)/len(ZP) - sum(ZP)/len(ZP)**2)
  ZP_mean = sum(weight*ZP)/sum(weight)
  ZP_wrms = np.sqrt(sum(weight*(ZP-ZP_mean)**2)/sum(weight))
  print 'Ensable APASS Uncertainties: ZP_rms, m_I_unc, m_APASS_unc',ZP_rms, np.sqrt(np.sum(np.array(m_I_unc)**2/len(m_I))), np.sqrt(np.sum(np.array(m_APASS_unc)**2/len(m_APASS)))
  #print ZP_mean, ZP_rms
  # WEIGHTED RMS
  # ESTIMATE THE WEIGHTED MEAN OF SEEING PARAMETERS
  beta_mean  = sum(weight*beta)/sum(weight)
  beta_wrms  = np.sqrt(sum(weight*(beta-beta_mean)**2)/sum(weight))

  # USE THE MEDIAN AND MEDIAN ABSOLUTE DEVIATION
  ZP_mean = np.median(ZP)
  ZP_wrms = 1.4826*np.median(np.abs(ZP - ZP_mean))
  ZP_rms =  1.4826*np.median(np.abs(ZP - ZP_mean))
  beta_mean  = np.median(beta)
  beta_unc  = 1.4826*np.median(np.abs(beta -  beta_mean))
  if(display>0):
  	# SORT BY INCREASING MAGNITUDE FOR EASE OF VIEWING
  	m_APASS_ord,m_I_ord = zip(*sorted(zip(m_APASS, m_I)))
  	m_APASS_ord,m_APASS_unc_ord = zip(*sorted(zip(m_APASS, m_APASS_unc)))
  	m_APASS_ord,m_I_unc_ord = zip(*sorted(zip(m_APASS, m_I_unc)))
  	m_APASS_ord,beta_ord = zip(*sorted(zip(m_APASS, beta)))
  	m_APASS_ord,beta_error_ord = zip(*sorted(zip(m_APASS, beta_error)))
  	m_APASS_ord,bkg_ord = zip(*sorted(zip(m_APASS, bkg)))
  	m_APASS_ord,bkg_error_ord = zip(*sorted(zip(m_APASS, bkg_error)))
	col='r'
	plt.figure(figsize=(8,12))
	if(filt=='Sloan_r'):
		plt.subplot(411)
		plt.errorbar(m_APASS_ord,np.array(m_APASS_ord)-np.array(m_I_ord),xerr=m_APASS_unc_ord, yerr=np.sqrt(np.array(m_I_unc_ord)**2+np.array(m_APASS_unc_ord)**2) ,fmt='.-', color=col)
		plt.title('Red Filter ZP')
 		plt.ylabel('Zero Point')
		#plt.plot()
	if(filt=='Sloan_g'):
		col='g'
		plt.subplot(411)
		plt.errorbar(m_APASS_ord,np.array(m_APASS_ord)-np.array(m_I_ord), xerr=m_APASS_unc_ord, yerr=np.sqrt(np.array(m_I_unc_ord)**2+np.array(m_APASS_unc_ord)**2),fmt='.-', color=col)
 		plt.ylabel('Zero Point')
 		plt.title('Green Filter ZP')
	plt.xlim(np.min(m_APASS_ord)-0.5, np.max(m_APASS_ord)+0.5)
	plt.grid(True)
	plt.subplot(412)
	plt.errorbar(m_APASS_ord, beta_ord, yerr=beta_error_ord, fmt='.-', color=col)
	#plt.plot(m_APASS_ord, beta_ord, '.-', color=col)
	plt.ylabel('beta')
	plt.title('Moffat beta')
	plt.grid(True)
	plt.xlim(np.min(m_APASS_ord)-0.5, np.max(m_APASS_ord)+0.5)
	plt.subplot(413)
	plt.errorbar(m_APASS_ord, bkg_ord, yerr=bkg_error_ord, fmt='.-', color=col)
	#plt.plot(m_APASS_ord, bkg_ord, '.-', color=col)
	plt.title('Fitted Nbkg')
	plt.ylabel('Counts')
	plt.grid(True)
	plt.subplots_adjust(left=0.2, right=0.95, hspace=0.4)
	plt.xlabel('APASS Source Magnitude', fontsize=14)
	plt.xlim(np.min(m_APASS_ord)-0.5, np.max(m_APASS_ord)+0.5)
	plt.suptitle(FM.fits_file_name.split('/')[-1], fontsize=18)
	plt.subplot(414)
	normed_parms = []
	for k in range(0,len(psf_parms)):
		normed_parms.append( psf_parms[k]/magnitude2flux(np.array(m_I[k]))*2.*np.sqrt(np.pi)*np.array(beta[k]) )
	med_parms = np.median(normed_parms, axis=0)
	plt.semilogy(np.abs(med_parms), 'k.-', alpha=0.3, lw=10 )    
	for k in range(0,len(psf_parms)):
		plt.semilogy(np.abs(normed_parms[k]), '.-')
	plt.savefig(out+'.png', dpi=50)
  '''

  # Estimate the median of the parameters
  normed_parms = []
  normed_parm_err = []
  for k in range(0,len(psf_parms)):
	#normed_parms.append( psf_parms[k]/magnitude2flux(np.array(m_I[k]))*2.*np.sqrt(np.pi)*np.array(beta[k]) )
	normed_parms.append( psf_parms[k]/np.array(S_CCD_list[k])*2.*np.sqrt(np.pi)*np.array(beta[k]) )
	normed_parm_err.append( psf_parm_err[k]/np.array(S_CCD_list[k])*2.*np.sqrt(np.pi)*np.array(beta[k]) )
  # take the median parameter value only for stars that pass the quality cut.
  normed_parms     = np.array(normed_parms)
  pass_quality_cut = np.array(pass_quality_cut)
  print 'pass_quality_cut', pass_quality_cut
  if np.sum(pass_quality_cut)==0:
    print 'NO SUCCESSFUL FITS TO PSF STARS'
    print 'EXITING PROGRAM'
    exit()
  med_parms = np.median(normed_parms[pass_quality_cut], axis=0)

  if(display>0):
	plt.figure()
	ax = plt.subplot(111)
	#ax.set_yscale('log')
	ax.set_xscale('log')
	plt.plot(range(1,len(normed_parms[k])+1),np.abs(med_parms), 'k.-', alpha=0.3, lw=10 )   
	mx = np.max(np.abs(normed_parms))*1.1
	for k in range(0,len(psf_parms)):
		#plt.plot(np.abs(normed_parms[k]), '.-')
		plt.errorbar(range(1,len(normed_parms[k])+1), np.abs(normed_parms[k]), yerr=normed_parm_err[k], fmt='.-')
        #print psf_parm_err[k]/normed_parms[k]
	nm_list = psf.get_nm(obj.nmax, obj.mmax)
	for i in range(0, len(nm_list)):
		n,m,imre = nm_list[i]
		if(n%2==0 and m==0):
			plt.plot([i+1,i+1],[-mx,mx], 'r--', lw=1) 
	plt.ylim(-mx, mx)
	plt.xlim(0.9, len(normed_parms[0])+0.5)
	plt.xlabel('Shapelet Index')
	plt.ylabel('Normalized Value')
	plt.title('APASS Star Normalized Shapelet Coefficients')
	plt.grid(True)
	plt.savefig(out+'_shapelet_coeffs.png', dpi=75)


  # HERE WE SHOULD ADD A COMBINED FIT
  # THE FIRST THING WOULD BE TO MASK PIXELS WITH >4 SIGMA RESIDUALS
  # INPUT THE FITTED BACKGROUND, POSITIONS, AND S_CCD FOR EACH STAR COMBINED WITH THE MEDIAN VALUE OF BETA AND THE SHAPELET COEFFICIENTS
  # FOR NOW, WE WILL JUST RETURN THE AVERAGE SHAPELET COEFFICIENTS
  beta_med  = np.median(beta)
  beta_unc  = 1.4826*np.median(np.abs(beta -  beta_med))
  print '\t\tBeta'.ljust(15), '%1.3f +/- %1.3f'%(beta_med,  beta_unc)
  return APASS_index_list, beta_med, beta_unc, med_parms, S_CCD_list, sig_S_CCD_list, chi_sq_list, max_chi_list

def readStarList(fnm):
	ra = []
	dec = []
	for line in file(fnm):
                if('#' not in line):
			        ra.append(float(line.split()[0]))
			        dec.append(float(line.split()[1]))
	return np.array(ra), np.array(dec)

def readStarList2(fnm):
	ra = []
	dec = []
	lc = 0
	data = np.genfromtxt(fnm)
	print data.shape
	for line in file(fnm):
	    print line
	    if('#' not in line and lc>0):
		    print line
		    ra.append(float(line.split(',')[0]))
		    dec.append(float(line.split(',')[1]))
	    lc+=1
	return np.array(ra), np.array(dec)

def starFit(FM, ra_star_list, dec_star_list,  beta, shapelet_coefficients, readnoise, nmax, mmax, N_px=31, display=0, outputFileTag='out_star'):
  #print '\nstarFit'

  nm_list = np.array(psf.get_nm(nmax, mmax))
  fn_indices = fn0_indices = np.where(nm_list[:,1]==0)[0]
  #print 'PSF_parms', len(PSF_parms)
  #print 'nm_list', len(nm_list)
  S_CCD = np.sum(shapelet_coefficients[fn0_indices])*(2*np.sqrt(np.pi)*beta)
  print '\tnormed_parameter S_CCD', S_CCD

  def shapeletMagnitude((x, y), x0, y0, bkg, sum_CCD):
    #p = sum_CCD*shapelet_coefficients
    #print '\t\t check S_CCD', sum_CCD, np.sum(p[fn0_indices])*(2*np.sqrt(np.pi)*beta)
    parameters = np.concatenate([[x0,y0, beta, bkg], sum_CCD*shapelet_coefficients])
    m = psf.shapelet_kernel(parameters, obj.nmax, obj.mmax, nx = len(x), ny = len(y))
    bkg = parameters[3]
    m = np.rot90(m)
    m = np.flipud(m)
    #m = amplitude*PSF_model + offset
    if(parameters[2]<0. or parameters[3]<0.): 
        m+=1.e9
    return m.ravel() + bkg

  x2d, y2d = np.mgrid[:N_px, :N_px]
  #star_index_list = []
  #chi_sq_list=[]
  #max_chi_list = []
  x0_list = []
  y0_list = []
  bkg_list = []
  S_CCD_list = []
  S_CCD_unc_list = []

  fit_success      = np.zeros(len(ra_star_list), dtype=np.bool)
  chi_sq_list      = np.zeros(len(ra_star_list))
  max_chi_list     = np.zeros(len(ra_star_list))
  x0_list          = np.zeros(len(ra_star_list))
  y0_list          = np.zeros(len(ra_star_list))
  bkg_list         = np.zeros(len(ra_star_list))
  S_CCD_list       = np.zeros(len(ra_star_list))
  S_CCD_unc_list   = np.zeros(len(ra_star_list))

  for k in range(0,len(ra_star_list)):
    print '\tStar Index =',k
    print '\t\tra, dec'.ljust(12), '= %1.5f , %1.5f'%(ra_star_list[k], dec_star_list[k])
    obj = SourceImage(FM, ra_star_list[k], dec_star_list[k], nmax, mmax, N_px)
    if(obj.image.shape != (N_px, N_px)): 
        print 'IMAGE FILE' 
        continue
    # get peak value and re-center around there
    x_max, y_max = np.unravel_index(obj.image.argmax(), obj.image.shape)
    #print 'x_max, y_max', x_max, y_max
    #print 'obj.image.shape', obj.image.shape
    # get the pixel coordinates for ra/dec,
    x,y = FM.bigw.wcs_world2pix(ra_star_list[k], dec_star_list[k],1)
    # add x_max, y_max to it
    x += y_max - N_px//2
    y += x_max - N_px//2
    # convert back to ra/dec
    ra, dec = FM.bigw.wcs_pix2world(x,y,1)
    # redefine obj in star-centered coordinates
    obj = SourceImage(FM, ra, dec, nmax, mmax, N_px)
    #print 'obj.image.shape', obj.image.shape
    if(obj.image.shape != (N_px, N_px)): continue
    x_max, y_max = np.unravel_index(obj.image.argmax(), obj.image.shape)
    #print x_max, y_max
    # guess values of x0, y0, and bkg
    x0_guess = float(x_max - (N_px- 1) / 2) 
    y0_guess = float(y_max - (N_px- 1) / 2)
    x2d, y2d = np.mgrid[:len(obj.image), :len(obj.image)]
    #print 'x_max, y_max', x_max, y_max
    #print 'x0_guess, y0_guess', x0_guess, y0_guess

    bkg_guess = np.median(obj.image)
    sum_CCD_guess = (np.max(obj.image)-bkg_guess)/1.266 # /(2.*np.sqrt(np.pi)*beta*1.226)
    guess_img = shapeletMagnitude( (x2d,y2d), x0_guess, y0_guess, bkg_guess, sum_CCD_guess)
    gain = FM.hdulist[0].header['GAIN']

    sigmas = np.sqrt( readnoise**2 + obj.image/gain )

    try:
      popt, pcov = opt.curve_fit(shapeletMagnitude, (x2d,y2d), obj.image.ravel(), sigma = sigmas.ravel(), absolute_sigma=True, p0=[x0_guess, y0_guess, bkg_guess, sum_CCD_guess])
    except:
      print '\t\tSTAR FIT FAILED'
      continue
    
    print '\t\tSTAR FIT SUCCEEDED'
    #print 'popt', popt
    #print 'pcov', pcov
    fitted_image = shapeletMagnitude( (x2d,y2d), *popt)
    df = obj.image - fitted_image.reshape(N_px,N_px)
    chi = df/sigmas
    print '\t\tS_CCD'.ljust(12), '= %1.0f +/- %1.0f'%( popt[3],np.sqrt(pcov[3][3]))
    print '\t\tchisq/ndof'.ljust(12), '= %1.3f'%(np.sum(chi**2)/float(chi.size-4))
    print '\t\tmax|chi|'.ljust(12), '= %1.3f'%np.max(np.abs(chi))

    fit_success[k]    = True
    chi_sq_list[k]    = np.sum(chi**2)/float(chi.size-4)
    max_chi_list[k]   = np.max(np.abs(chi))
    x0_list[k]        = popt[0]
    y0_list[k]        = popt[1]
    bkg_list[k]       = popt[2]
    S_CCD_list[k]     = popt[3]
    S_CCD_unc_list[k] = np.sqrt(pcov[3][3])
    

    #star_index_list.append(k)
    #chi_sq_list.append(np.sum(chi**2)/float(chi.size-4))
    #max_chi_list.append(np.max(np.abs(chi)))
    #S_CCD_list.append(popt[3])
    #S_CCD_unc_list.append(np.sqrt(pcov[3][3]))
    if(display>2):
  	  plt.figure()
  	  mx = np.max([np.max(obj.image), np.max(fitted_image)])
  	  mn = np.min([np.min(obj.image), np.min(fitted_image)])
  	  plt.subplot(222)
  	  plt.imshow(obj.image, interpolation='none', vmin=mn, vmax=mx, cmap='viridis')
  	  plt.colorbar()
  	  plt.title('data')
  	  #plt.subplot(232)
  	  #plt.imshow(guess_img.reshape(N_px,N_px), interpolation='none', cmap='viridis')
  	  #plt.colorbar()
  	  #plt.subplot(233)
  	  #plt.imshow(obj.image-guess_img.reshape(N_px,N_px), interpolation='none', cmap='viridis')
  	  #plt.colorbar()
  	  plt.subplot(221)
  	  plt.imshow(fitted_image.reshape(N_px,N_px), interpolation='none', vmin=mn, vmax=mx, cmap='viridis')
  	  plt.title('fit')
  	  plt.colorbar()
  	  plt.subplot(223)
  	  plt.imshow(df, interpolation='none', vmin=-np.max(np.abs(df)), vmax=np.max(np.abs(df)), cmap='seismic')
  	  plt.colorbar()
  	  plt.title('residual')
  	  plt.subplot(224)
  	  plt.imshow(chi, interpolation='none', vmin = -np.max(np.abs(chi)), vmax = np.max(np.abs(chi)), cmap='seismic')
  	  plt.colorbar()
  	  plt.title('chi')
  	  plt.savefig(outputFileTag+'_star_%d.png'%k)
  return fit_success, chi_sq_list, max_chi_list, x0_list, y0_list, bkg_list, S_CCD_list, S_CCD_unc_list



def quadFit(FM, ra_qsr, dec_qsr, ra_images, dec_images, ra_lens, dec_lens, readnoise, beta, shapelet_coeffs, nmax, mmax, N_px, galFit=False, gal_S_CCD = None, outputFileTag='out', display=0, emcee_level=1):
  print '\nquadFit'
  # GET THE QUASAR IMAGE
  obj = SourceImage(FM, ra_qsr, dec_qsr, nmax, mmax, N_px)

  N_bkg = np.median(obj.image)
  #intg, intg_unc = estimate_total_light(obj,N_bkg, FM.readnoise, display=False)

  # DETERMINE THE LEFT RIGHT ORIENTATION OF THE IMAGE BASED ON ITS RIGHT ASCENCION DIFFERENCE WITH PIXEL VALUES
  xv,yv = FM.bigw.wcs_world2pix(ra_qsr,dec_qsr,1)
  r0,d0 = FM.bigw.wcs_pix2world(xv,yv,1)
  r1,d1 = FM.bigw.wcs_pix2world(xv+1,yv+1,1)
  fl = True
  if(r1-r0<0): fl =False

  # SET AN INITIAL AMPLITUDE SCALE, THESE VALUES ARE APPROXIMATED FROM Kochenek, 2006
  amp_scale = 1. # This value will be re-normalized after correlation.
  x0,y0 = 0., 0. 
  amp0 = amp_scale    
  amp1 = amp_scale/1.2
  amp2 = amp_scale/1.5
  amp3 = amp_scale/2.
  amp_lens = amp_scale*0. # zero this out in the initial fit

  # convert image ra and dec from arcseconds to degrees
  dec_images = dec_qsr + dec_images/3600.
  ra_images  = ra_qsr  + ra_images/3600.
  dec_lens   = dec_qsr + dec_lens/3600.
  ra_lens    = ra_qsr  + ra_lens/3600.

  # get the image positions in pixels
  x_images,y_images = FM.bigw.wcs_world2pix(ra_images,dec_images,1)
  x_images -= xv #- (N_px-1)/2
  y_images -= yv #- (N_px-1)/2

  # set lower and upper bounds to the image pixel positions 
  # multiplied by 1.e-5 because I don't want these to vary
  x_imagesU = x_images - 1.e-5*3.*0.003/float(FM.hdulist[0].header['PIXSCALE']) # 3-sigma bound
  x_imagesL = x_images + 1.e-5*3.*0.003/float(FM.hdulist[0].header['PIXSCALE']) # 3-sigma bound
  y_imagesU = y_images - 1.e-5*3.*0.003/float(FM.hdulist[0].header['PIXSCALE']) # 3-sigma bound
  y_imagesL = y_images + 1.e-5*3.*0.003/float(FM.hdulist[0].header['PIXSCALE']) # 3-sigma bound
  for k in range(0,len(x_imagesU)):
    if x_imagesU[k]<x_imagesL[k]:
        x_imagesU[k], x_imagesL[k] = x_imagesL[k], x_imagesU[k]
    if y_imagesU[k]<y_imagesL[k]:
        y_imagesU[k], y_imagesL[k] = y_imagesL[k], y_imagesU[k]

  x_lens,y_lens = FM.bigw.wcs_world2pix(ra_lens,dec_lens,1)
  # x_lens and y_lens are offsets from the lens coordinates
  x_lens -= xv #- (N_px-1)/2
  y_lens -= yv #- (N_px-1)/2

  x_lensL = x_lens - 1.e-5*3.*0.003/float(FM.hdulist[0].header['PIXSCALE']) # 3-sigma bound
  x_lensU = x_lens + 1.e-5*3.*0.003/float(FM.hdulist[0].header['PIXSCALE']) # 3-sigma bound
  y_lensL = y_lens - 1.e-5*3.*0.003/float(FM.hdulist[0].header['PIXSCALE']) # 3-sigma bound 
  y_lensU = y_lens + 1.e-5*3.*0.003/float(FM.hdulist[0].header['PIXSCALE']) # 3-sigma bound
  y_lensL = y_lens - 1.e-5*3.*0.003/float(FM.hdulist[0].header['PIXSCALE']) # 3-sigma bound 

  # Check the ordering of upper and lower. 
  if(x_lensL>x_lensU):
    x_lensL, x_lensU = x_lensU, x_lensL
  if(y_lensL>y_lensU):
    y_lensL, y_lensU = y_lensU, y_lensL

  # lens deVaucouleurs radius
  r0_lens  = 1.5 # arcseconds, from Courbin+ arXiv:1009.1473
  r0_lensU = 1.5+3.*0.08 # 3-sigma upped bound from Courbin+ arXiv:1009.1473
  r0_lensL = 1.5-3.*0.08 # 3-sigma upped bound from Courbin+ arXiv:1009.1473
  # convert to pixel-space values
  r0_lens  /= float(FM.hdulist[0].header['PIXSCALE'])
  r0_lensU /= float(FM.hdulist[0].header['PIXSCALE'])
  r0_lensL /= float(FM.hdulist[0].header['PIXSCALE'])

  # lens ellipticity, from Courbin+ arXiv:1009.1473
  el_lens = 0.3 # 0.09
  err_el =  0.01     # 0.01
  q_lens  = 1./(1.-el_lens)
  q_lensU = q_lens + 3.*err_el/(1.-el_lens)**2 # 3-sigma upped bound from Courbin+ arXiv:1009.1473
  q_lensL = q_lens - 3.*err_el/(1.-el_lens)**2 # 3-sigma upped bound from Courbin+ arXiv:1009.1473

  posang_lens  = 174.8
  posang_lensU = posang_lens + 3.*1.7 # 3-sigma upped bound from Courbin+ arXiv:1009.1473
  posang_lensL = posang_lens - 3.*1.7 # 3-sigma upped bound from Courbin+ arXiv:1009.1473
  #print xv, x_images, x_lens
  #print yv, y_images, y_lens

  # create x,y grid
  xg, yg = np.mgrid[:N_px,:N_px] 
  # DEFINE ARRAY OF INPUT PARAMETERS
  # PRODUCE A MODELED IMAGE
  #qim = obj.quad_image_model(x0,y0,x_images,y_images, amp0,amp1,amp2,amp3, alpha, beta, N_bkg, len(obj.image), fl, x_lens, y_lens, amp_lens, r0_lens, q_lens, posang_lens)

  # FOR THE CORRELATION AND FIRST IMAGE ALIGNMENT, WE WON'T INCLUDE THE LENSING GALAXY
  # set the background to 0 for cross-correlation
  qim = obj.quad_image_model((xg,yg, beta, shapelet_coeffs, len(obj.image), fl),x0,y0,x_images[0],x_images[1],x_images[2],x_images[3],
                             y_images[0],y_images[1],y_images[2],y_images[3], 
                             amp0,amp1,amp2,amp3, 0.).reshape(N_px, N_px)
  '''
  plt.figure()
  plt.subplot(211)
  plt.imshow(qim, interpolation='none', cmap='YlGnBu_r')
  plt.colorbar()
  plt.subplot(212)
  plt.imshow(obj.image, interpolation='none', cmap='YlGnBu_r')
  plt.colorbar()
  plt.savefig('qim_test.png')
  #exit()
  '''
  # CROSS CORRELATE THE DATA TO THE MODEL TO FIND OUT WHERE IT IS ON THE PIXEL GRID
  corr = signal.correlate2d(obj.image, qim, boundary='symm', mode='same')
  # GET THE LOCATION OF THE CORRELATION PEAK
  corr_x_max, corr_y_max = np.unravel_index(corr.argmax(), corr.shape)
  # GET THE IMAGE LOCATION CORRECTIONS
  dx = corr_x_max-15
  dy = corr_y_max-15
  # SOME FUNINESS HERE. I DETERMINED IT FROM SCANNING THROUGH SEVERAL IMAGES
  #print 'corr_x_max, corr_y_max', corr_x_max, corr_y_max
  #print 'dx, dy', dx,dy
  #if(corr_x_max>=15): dx = corr_x_max - 31
  #if(corr_y_max>=15): dy = corr_y_max - 31
  # DEFINE NEW PARAMETERS FOR THE MODELED IMAGE
  x0 += dx
  y0 += dy
  amp_fac = np.max(obj.image - N_bkg)/np.max(qim)
  amp0 *= amp_fac
  amp1 *= amp_fac
  amp2 *= amp_fac
  amp3 *= amp_fac
  amp_lens *= amp_fac
  #print 'amp correction factor', 
  # PRODUCE THE MODEL IMAGE WITHOUT LENSING GALAXY
  qim2 = obj.quad_image_model((xg,yg, beta, shapelet_coeffs, len(obj.image), fl), x0,y0,x_images[0],x_images[1],x_images[2],x_images[3],
                             y_images[0],y_images[1],y_images[2],y_images[3], 
                             amp0,amp1,amp2,amp3, N_bkg ).reshape(N_px, N_px)
  '''
  print 'dx, dy', dx,dy
  plt.figure()
  plt.subplot(221)
  plt.imshow(qim2, interpolation='none', cmap='viridis')
  plt.title('model')
  plt.colorbar()
  plt.subplot(222)
  plt.imshow(obj.image, interpolation='none', cmap='viridis')
  plt.title('data')
  plt.colorbar()
  plt.subplot(223)
  plt.imshow(corr, interpolation='none', cmap='viridis')
  plt.title('correlation')
  plt.colorbar()
  plt.subplot(224)
  plt.imshow(qim2, interpolation='none', cmap='viridis')
  plt.title('model')
  plt.colorbar()
  plt.savefig('qim2_test.png')
  exit()
  '''
 
  #### FIT ONCE WITHOUT THE LENSING GALAXY ############################################################

  # curve_fit qim
  gain = FM.hdulist[0].header['GAIN']
  sigmas = np.sqrt( readnoise**2 + obj.image/gain )

  '''
  p0 = [x0,  y0,  x_images[0], x_images[1], x_images[2], x_images[3], y_images[0], y_images[1], y_images[2], y_images[3], amp0,   amp1,   amp2, amp3,  N_bkg]
  pU = [15., 15., x_imagesU[0],x_imagesU[1],x_imagesU[2],x_imagesU[3],y_imagesU[0],y_imagesU[1],y_imagesU[2],y_imagesU[3], np.inf,np.inf,np.inf,np.inf,np.inf]
  pL = [-15.,-15.,x_imagesL[0],x_imagesL[1],x_imagesL[2],x_imagesL[3],y_imagesL[0],y_imagesL[1],y_imagesL[2],y_imagesL[3], 0.,    0.,    0.,    0.,    0.,   ]
  bounds = [pL,pU]
  print '\tfitting without the lensing galaxy'
  popt_ng, pcov_ng = opt.curve_fit(obj.quad_image_model, (xg, yg, beta, shapelet_coeffs, len(obj.image), fl), obj.image.ravel(), sigma = sigmas.ravel(), absolute_sigma=True, p0=p0, bounds=bounds)
  '''

  p0 = [x0,    y0,   amp0,   amp1,   amp2,   amp3,   N_bkg]
  pU = [15.,   15.,  np.inf, np.inf, np.inf, np.inf, np.inf]
  pL = [-15., -15.,  0.,     0.,     0.,     0.,     0.,   ]
  bounds = [pL,pU]

  amp_lens=0.
  if(gal_S_CCD != None):
    amp_lens = gal_S_CCD
  # wrapper to fix the image positions relative to the field position
  def qmi_ng( (xg,yg, beta, shapelet_coeffs, N_pix, flip,  _x1, _x2, _x3, _x4, _y1, _y2, _y3, _y4, _xlens, _ylens, _amp_lens, _r0_lens, _q_lens, _posang_lens), x0, y0, amp0, amp1, amp2, amp3, N_bkg):
    return obj.quad_image_model((xg,yg, beta, shapelet_coeffs, N_pix, flip), x0, y0, _x1, _x2, _x3, _x4, _y1, _y2, _y3, _y4, amp0, amp1, amp2, amp3, N_bkg, _xlens, _ylens, _amp_lens, _r0_lens, _q_lens, _posang_lens)

  print '\tFitting w/o lensing galaxy'
  popt_ng, pcov_ng = opt.curve_fit(qmi_ng, (xg, yg, beta, shapelet_coeffs, len(obj.image), fl, x_images[0], x_images[1], x_images[2], x_images[3], y_images[0], y_images[1], y_images[2], y_images[3], x_lens, y_lens, amp_lens, r0_lens, q_lens, posang_lens), obj.image.ravel(), sigma = sigmas.ravel(), absolute_sigma=True, p0=p0, bounds=bounds)


  print '\tParameters'
  print '\t\tx0    = %+1.2e +/- %+1.2e'%(popt_ng[0], np.sqrt(pcov_ng[0,0]))
  print '\t\ty0    = %+1.2e +/- %+1.2e'%(popt_ng[1], np.sqrt(pcov_ng[1,1]))
  print '\t\tamp_0 = %+1.2e +/- %+1.2e'%(popt_ng[2], np.sqrt(pcov_ng[2,2]))
  print '\t\tamp_1 = %+1.2e +/- %+1.2e'%(popt_ng[3], np.sqrt(pcov_ng[3,3]))
  print '\t\tamp_2 = %+1.2e +/- %+1.2e'%(popt_ng[4], np.sqrt(pcov_ng[4,4]))
  print '\t\tamp_3 = %+1.2e +/- %+1.2e'%(popt_ng[5], np.sqrt(pcov_ng[5,5]))
  print '\t\tN_bkg = %+1.2e +/- %+1.2e'%(popt_ng[6], np.sqrt(pcov_ng[6,6]))
  #print popt_ng
  qim3 = qmi_ng((xg,yg, beta, shapelet_coeffs, len(obj.image), fl, x_images[0], x_images[1], x_images[2], x_images[3], y_images[0], y_images[1], y_images[2], y_images[3], x_lens, y_lens, amp_lens, r0_lens, q_lens, posang_lens), *popt_ng).reshape(N_px, N_px)

  # ESTIMATE ITS DISTRIBUTION OF CHI VALUES
  res = obj.image-qim3
  chi  = res/sigmas
  chi_x_max, chi_y_max = np.unravel_index(np.abs(chi).argmax(), chi.shape)
  chi2 = np.sum(chi**2)/(float(obj.image.size - len(popt_ng)))
  print '\t\tchi^2 = %+1.2f'%(chi2)
  print '\t\tm_chi = %+1.2f'%(np.max(np.abs(chi)))
  #print 'quadFit chi2, max_chi', chi2, np.max(np.abs(chi)), chi_x_max, chi_y_max

  #'''
  #print 'dx, dy', dx,dy
  mn = np.min(obj.image)
  mx = np.max(obj.image)
  if(display>2):
    plt.figure()
    plt.subplot(221)
    plt.imshow(qim3, interpolation='none', vmin = mn, vmax = mx, cmap='viridis')
    plt.title('fit')
    plt.colorbar()
    plt.subplot(222)
    plt.imshow(obj.image, interpolation='none', vmin = mn, vmax = mx, cmap='viridis')
    plt.title('data')
    plt.colorbar()
    plt.subplot(223)
    plt.imshow(res, interpolation='none', vmin=-np.max(np.abs(res)), vmax=np.max(np.abs(res)), cmap='seismic')
    plt.title('residuals')
    plt.colorbar()
    plt.subplot(224)
    plt.imshow(chi, interpolation='none', vmin=-np.max(np.abs(chi)), vmax=np.max(np.abs(chi)), cmap='seismic')
    plt.plot([chi_y_max],[chi_x_max], 'o', ms=8, mfc='none', mec='g', mew=2 )
    plt.xlim(0-0.5,N_px-0.5)
    plt.ylim(N_px-0.5, 0-0.5)
    plt.title('chi')
    plt.colorbar()
    #plt.suptitle()
    plt.savefig('%s_quadFit_NG.png'%outputFileTag)
  #'''
  if(galFit==False):
    print '\tgalFit',galFit
    print '\tamp_lens',amp_lens
    '''
    amp0_fit  = popt_ng[10]
    amp1_fit  = popt_ng[11]
    amp2_fit  = popt_ng[12]
    amp3_fit  = popt_ng[13]
    N_bkg_fit = popt_ng[14]
    amp0_err  = np.sqrt(pcov_ng[10][10])
    amp1_err  = np.sqrt(pcov_ng[11][11])
    amp2_err  = np.sqrt(pcov_ng[12][12])
    amp3_err  = np.sqrt(pcov_ng[13][13])
    N_bkg_err = np.sqrt(pcov_ng[14][14])
    '''
    #print pcov_ng[10:][10:]
    return popt_ng, pcov_ng, chi2, np.max(np.abs(chi))


  # print 'Fitting the central galaxy'
  ##################################################################################################

  # DO NOT LET LENS PARAMETERS VARY, EXCEPT FOR AMPLITUDE
  # CREATE A WRAPPED WHEN THE LENS PARAMETERS ARE HELD FIXED
  # IF A ZERO POINT IS PASSED, WE WANT TO BOUND THE AMPLITUDE.

  # start with post-fit values
  x0,  y0, amp0,  amp1,  amp2,  amp3,  N_bkg = popt_ng

  def qim_wg( (xg,yg, beta, shapelet_coeffs, N_pix, flip,  _x1, _x2, _x3, _x4, _y1, _y2, _y3, _y4, _xlens, _ylens, _r0_lens, _q_lens, _posang_lens), x0, y0, amp0, amp1, amp2, amp3, N_bkg, amp_lens):
    return obj.quad_image_model((xg,yg, beta, shapelet_coeffs, N_pix, flip), x0, y0, _x1, _x2, _x3, _x4, _y1, _y2, _y3, _y4, amp0, amp1, amp2, amp3, N_bkg, _xlens, _ylens, amp_lens, _r0_lens, _q_lens, _posang_lens)
  # SET THE INITIAL LENS AMPLITUDE TO 1/10th THE LEVEL OF THE BACKGROUND
  print ''
  print '\tLensing Galaxy Parameters'
  print '\t\tr0_lens %1.2f'%r0_lens
  print '\t\tq_lens  %1.2f'%q_lens
  print '\t\tposang  %1.2f'%posang_lens
  amp_lens0 = 1.e9
  m = deVaucouleurs_model(x0+(N_px-1)/2, y0+(N_px-1)/2, amp_lens0, r0_lens, beta, shapelet_coeffs, nmax, mmax, q = q_lens, posang=posang_lens, npx=31)
  amp_lens0 *= 0.1*N_bkg/np.max(m)
  m = deVaucouleurs_model(x0+(N_px-1)/2, y0+(N_px-1)/2, amp_lens0, r0_lens, beta, shapelet_coeffs, nmax, mmax, q = q_lens, posang=posang_lens, npx=31)
  amp_lens0U = np.inf
  amp_lens0L = 0.

  # perform combined fit
  '''
  p0 = [x0,  y0,  x_images[0], x_images[1], x_images[2], x_images[3], y_images[0], y_images[1], y_images[2], y_images[3],  amp0,  amp1,  amp2,  amp3,  N_bkg,  x_lens, y_lens, amp_lens0, r0_lens, q_lens, posang_lens]
  pU = [15., 15., x_imagesU[0],x_imagesU[1],x_imagesU[2],x_imagesU[3],y_imagesU[0],y_imagesU[1],y_imagesU[2],y_imagesU[3], np.inf,np.inf,np.inf,np.inf,np.inf, x_lensU,y_lensU,amp_lens0U,r0_lensU,q_lensU,posang_lensU]
  pL = [-15.,-15.,x_imagesL[0],x_imagesL[1],x_imagesL[2],x_imagesL[3],y_imagesL[0],y_imagesL[1],y_imagesL[2],y_imagesL[3], 0.,    0.,    0.,    0.,    0.,     x_lensL,y_lensL,amp_lens0L,r0_lensL,q_lensL,posang_lensL]
  print '\tfitting with the lensing galaxy'
  bounds = [pL,pU]
  popt_wg, pcov_wg = opt.curve_fit(obj.quad_image_model, (xg, yg, beta, shapelet_coeffs, len(obj.image), fl), obj.image.ravel(), sigma = sigmas.ravel(), absolute_sigma=True, p0=p0, bounds=bounds)
  '''
  p0 = [x0,   y0,  amp0,   amp1,   amp2,   amp3,   N_bkg,  amp_lens0]
  pU = [15.,  15., np.inf, np.inf, np.inf, np.inf, np.inf, amp_lens0U]
  pL = [-15.,-15., 0.,     0.,     0.,     0.,     0.,     amp_lens0L]
  print ''
  print '\tFitting w/ lensing galaxy'
  bounds = [pL,pU]
  popt_wg, pcov_wg = opt.curve_fit(qim_wg, (xg, yg, beta, shapelet_coeffs, len(obj.image), fl, x_images[0], x_images[1], x_images[2], x_images[3], y_images[0], y_images[1], y_images[2], y_images[3], x_lens, y_lens, r0_lens, q_lens, posang_lens), obj.image.ravel(), sigma = sigmas.ravel(), absolute_sigma=True, p0=p0, bounds=bounds)

  print '\tParameters'
  print '\t\tx0    = %+1.2e +/- %+1.2e'%(popt_wg[0], np.sqrt(pcov_wg[0,0]))
  print '\t\ty0    = %+1.2e +/- %+1.2e'%(popt_wg[1], np.sqrt(pcov_wg[1,1]))
  print '\t\tamp_0 = %+1.2e +/- %+1.2e'%(popt_wg[2], np.sqrt(pcov_wg[2,2]))
  print '\t\tamp_1 = %+1.2e +/- %+1.2e'%(popt_wg[3], np.sqrt(pcov_wg[3,3]))
  print '\t\tamp_2 = %+1.2e +/- %+1.2e'%(popt_wg[4], np.sqrt(pcov_wg[4,4]))
  print '\t\tamp_3 = %+1.2e +/- %+1.2e'%(popt_wg[5], np.sqrt(pcov_wg[5,5]))
  print '\t\tN_bkg = %+1.2e +/- %+1.2e'%(popt_wg[6], np.sqrt(pcov_wg[6,6]))
  print '\t\tamp_L = %+1.2e +/- %+1.2e'%(popt_wg[7], np.sqrt(pcov_wg[7,7]))

  # produce best fit models
  qim_gal = qim_wg((xg,yg, beta, shapelet_coeffs, len(obj.image), fl, x_images[0], x_images[1], x_images[2], x_images[3], y_images[0], y_images[1], y_images[2], y_images[3], x_lens, y_lens, r0_lens, q_lens, posang_lens), *popt_wg).reshape(N_px, N_px)
  p_gal = popt_wg.copy()
  p_gal[2:7] *= 0.
  g = qim_wg((xg,yg, beta, shapelet_coeffs, len(obj.image), fl, x_images[0], x_images[1], x_images[2], x_images[3], y_images[0], y_images[1], y_images[2], y_images[3], x_lens, y_lens, r0_lens, q_lens, posang_lens), *p_gal).reshape(N_px, N_px)

  # calculate residuals, chi, and chi squared values
  res = obj.image-qim_gal
  chi  = res/sigmas
  chi_x_max, chi_y_max = np.unravel_index(np.abs(chi).argmax(), chi.shape)
  chi2 = np.sum(chi**2)/(float(obj.image.size - len(popt_ng)))

  print '\t\tchi^2 = %+1.2f'%(chi2)
  print '\t\tm_chi = %+1.2f'%(np.max(np.abs(chi)))

  #print 'dx, dy', dx,dy
  if(display>2):
    mn = np.min(obj.image)
    mx = np.max(obj.image)
    plt.figure()
    plt.subplot(221)
    plt.imshow(qim_gal, interpolation='none', vmin = mn, vmax = mx, cmap='viridis')
    plt.title('fit')
    plt.colorbar()
    plt.subplot(222)
    plt.imshow(obj.image, interpolation='none', vmin = mn, vmax = mx, cmap='viridis')
    plt.title('data')
    plt.colorbar()
    plt.subplot(223)
    plt.imshow(res, interpolation='none', vmin=-np.max(np.abs(res)), vmax=np.max(np.abs(res)), cmap='seismic')
    plt.title('residuals')
    plt.colorbar()
    plt.subplot(224)
    plt.imshow(chi, interpolation='none', vmin=-np.max(np.abs(chi)), vmax=np.max(np.abs(chi)), cmap='seismic')
    plt.plot([chi_y_max],[chi_x_max], 'o', ms=8, mfc='none', mec='g', mew=2 )
    plt.xlim(0-0.5,N_px-0.5)
    plt.ylim(N_px-0.5, 0-0.5)
    plt.title('chi')
    plt.colorbar()
    #plt.suptitle()
    plt.savefig('%s_quadFit_WG.png'%outputFileTag)

    '''
    plt.figure()
    plt.subplot(111)
    plt.imshow(g, interpolation='none', cmap='viridis')
    plt.colorbar()
    plt.contour(obj.image, levels=np.linspace(0.5*np.max(obj.image), np.max(obj.image), 10))
    plt.savefig('%s_fitgal.png'%outputFileTag)

    plt.figure()
    plt.subplot(111)
    plt.imshow(g, interpolation='none', cmap='viridis')
    plt.colorbar()
    plt.contour(qim_gal, levels=np.linspace(0.5*np.max(qim_gal), np.max(qim_gal), 10))
    plt.savefig('%s_fitgal2.png'%outputFileTag)
    '''

    plt.figure()
    plt.subplot(111)
    plt.imshow(obj.image, interpolation='none', cmap='viridis')
    plt.colorbar()
    plt.contour(g, levels=np.linspace(0.5*np.max(g), 0.51*np.max(g), 1), colors='r')
    p_im1 = popt_wg.copy()
    #p_im1[[11,12,13,14,17]] *= 0.
    #im1 = obj.quad_image_model((xg,yg, beta, shapelet_coeffs, len(obj.image), fl), *p_im1).reshape(N_px, N_px)
    p_im1[[3,4,5,6,7]] *= 0.
    im1 = qim_wg((xg,yg, beta, shapelet_coeffs, len(obj.image), fl, x_images[0], x_images[1], x_images[2], x_images[3], y_images[0], y_images[1], y_images[2], y_images[3], x_lens, y_lens, r0_lens, q_lens, posang_lens), *p_im1).reshape(N_px, N_px)
    plt.contour(im1, levels=np.linspace(0.5*np.max(im1), 0.51*np.max(im1), 1), colors='r')

    p_im2 = popt_wg.copy()
    p_im2[[2,4,5,6,7]] *= 0.
    im2 = qim_wg((xg,yg, beta, shapelet_coeffs, len(obj.image), fl, x_images[0], x_images[1], x_images[2], x_images[3], y_images[0], y_images[1], y_images[2], y_images[3], x_lens, y_lens, r0_lens, q_lens, posang_lens), *p_im2).reshape(N_px, N_px)
    plt.contour(im2, levels=np.linspace(0.5*np.max(im2), 0.51*np.max(im2), 1), colors='r')

    p_im3 = popt_wg.copy()
    p_im3[[2,3,5,6,7]] *= 0.
    im3 = qim_wg((xg,yg, beta, shapelet_coeffs, len(obj.image), fl, x_images[0], x_images[1], x_images[2], x_images[3], y_images[0], y_images[1], y_images[2], y_images[3], x_lens, y_lens, r0_lens, q_lens, posang_lens), *p_im3).reshape(N_px, N_px)
    plt.contour(im3, levels=np.linspace(0.5*np.max(im3), 0.51*np.max(im3), 1), colors='r')

    p_im4 = popt_wg.copy()
    p_im4[[2,3,4,6,7]] *= 0.
    im4 = qim_wg((xg,yg, beta, shapelet_coeffs, len(obj.image), fl, x_images[0], x_images[1], x_images[2], x_images[3], y_images[0], y_images[1], y_images[2], y_images[3], x_lens, y_lens, r0_lens, q_lens, posang_lens), *p_im4).reshape(N_px, N_px)
    plt.contour(im4, levels=np.linspace(0.5*np.max(im4), 0.51*np.max(im4), 1), colors='r')
    plt.savefig('%s_fitgal_contours.png'%outputFileTag)

  return popt_wg, pcov_wg, chi2, np.max(np.abs(chi))

# THE FUNCTION BELOW NEEDS MAJOR UPDATES! SEE quadFit ABOVE.
def emceeQuadFit(FM, ra_qsr, dec_qsr, ZP_mean, ZP_rms, alpha, beta, alpha_err, beta_err, alpha_beta_corr_coeff, m1, m2, m3, m4, N_px, ra_images, dec_images, outputFileTag='out'):
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

  seeing_cov = np.array([[alpha_err**2, alpha_beta_corr_coeff*alpha_err*beta_err],
			 [ alpha_beta_corr_coeff*alpha_err*beta_err, beta_err**2]])
  det = np.linalg.det(seeing_cov)
  seeing_icov = np.linalg.inv(seeing_cov)
  print 'seeing_cov', seeing_cov
  print 'seeing_det', det
  print 'seeing_icov', seeing_icov
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
          seeing_diff = np.array([alpha-_alpha, beta-_beta])
          lnp = -np.dot( seeing_diff, np.dot(seeing_icov,seeing_diff) )/2.0
	  return lnp
    else:
        return -np.inf

  def lnprob(theta, obj, N_pix, flip):
    lp = lnprior(theta)
    #print lp, theta
    if not np.isfinite(lp):
        return -np.inf
    mll = obj.moffat_chi_sq(theta, N_pix, flip, x_images, y_images)
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
  # convert image ra and dec from arcseconds to degrees
  dec_images = dec_qsr + dec_images/3600.
  ra_images = ra_qsr   + ra_images/3600.
  print dec_images
  print ra_images
  x_images,y_images = FM.bigw.wcs_world2pix(ra_images,dec_images,1)
  x_images -= xv - (N_px-1)/2
  y_images -= yv - (N_px-1)/2
  print xv, x_images
  print yv, y_images

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
  fig= corner.corner(samples, labels=labels)
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

  fig= corner.corner(samples, labels=labels)
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

##########################################

def emceeQuadFitLensGal(FM, ra_qsr, dec_qsr, ZP_mean, ZP_rms, alpha, beta, alpha_err, beta_err, alpha_beta_corr_coeff, m1, m2, m3, m4, N_px, ra_images, dec_images, ra_lens, dec_lens, outputFileTag='out'):
  print '\n\t####################################'
  print   '\t# STARTING EMCEE ###################'
  print   '\t####################################\n'
  obj = SourceImage(FM, ra_qsr, dec_qsr, N_px)
  N_bkg_guess = np.median(obj.image)

  x0,y0, amp1, amp2, amp3, amp4, amp_lens, alpha_qsr, beta_qsr, N_bkg_qsr = FM.qsr_min_parms

  # DETERMINE THE IMAGE ORIENTATION
  xv,yv = FM.bigw.wcs_world2pix(ra_qsr,dec_qsr,1)
  r0,d0 = FM.bigw.wcs_pix2world(xv,yv,1)
  r1,d1 = FM.bigw.wcs_pix2world(xv+1,yv+1,1)
  fl = True
  if(r1-r0<0): fl =False
  # convert image ra and dec from arcseconds to degrees
  dec_images = dec_qsr + dec_images/3600.
  ra_images = ra_qsr   + ra_images/3600.
  print dec_images
  print ra_images
  x_images,y_images = FM.bigw.wcs_world2pix(ra_images,dec_images,1)
  x_images -= xv - (N_px-1)/2
  y_images -= yv - (N_px-1)/2
  x_lens,y_lens = FM.bigw.wcs_world2pix(ra_lens,dec_lens,1)
  x_lens -= xv - (N_px-1)/2
  y_lens -= yv - (N_px-1)/2
  r0_lens = 1.5/float(FM.hdulist[0].header['PIXSCALE'])
  el_lens = 0.09
  q_lens = 1./np.sqrt(1.-el_lens**2)
  posang_lens = 174.8
  print xv, x_images, x_lens
  print yv, y_images, y_lens

  # INITIAL ESTIMATE OF LENS AMPLITUDE ASSUMING IT IS IN THE SCALE OF THE NOISE
  #print 'amp_lens initialize', N_bkg_qsr, FM.readnoise
  #amp_lens = 126669.015793*np.sqrt(N_bkg_qsr + FM.readnoise**2)*r0_lens**2


  #amp1 = 10**((m1-ZP_mean)/(-2.5))
  #amp2 = 10**((m2-ZP_mean)/(-2.5))
  #amp3 = 10**((m3-ZP_mean)/(-2.5))
  #amp4 = 10**((m4-ZP_mean)/(-2.5))
  #N_bkg = np.median(obj.image)

  #obj = SourceImage(FM, ra_qsr, dec_qsr, N_px)

  seeing_cov = np.array([[alpha_err**2, alpha_beta_corr_coeff*alpha_err*beta_err],
			 [ alpha_beta_corr_coeff*alpha_err*beta_err, beta_err**2]])
  det = np.linalg.det(seeing_cov)
  seeing_icov = np.linalg.inv(seeing_cov)
  print 'seeing_cov', seeing_cov
  print 'seeing_det', det
  print 'seeing_icov', seeing_icov
  def lnprior(_theta):
    _x0, _y0, \
    _amplitude0, \
    _amplitude1, \
    _amplitude2, \
    _amplitude3, \
    _amp_lens, \
    _r0_lens, \
    _q_lens, \
    _posang_lens, \
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
	and _r0_lens>1.35/float(FM.hdulist[0].header['PIXSCALE']) and _r0_lens<1.65/float(FM.hdulist[0].header['PIXSCALE']) \
	and _q_lens>1.002 and _q_lens<1.006 \
	and _posang_lens>171.4 and _posang_lens<178.2\
	and _alpha>0. and _alpha < 100.*alpha\
	and _beta>0. and  _beta  < 100.*beta \
	and _N_bkg>0. and _N_bkg < 100.*N_bkg_guess:
          seeing_diff = np.array([alpha-_alpha, beta-_beta])
          lnp = -np.dot( seeing_diff, np.dot(seeing_icov,seeing_diff) )/2.0
          #print 'lnp', lnp
	  return lnp
    else:
        '''
        print '\t',np.abs(_x0),N_px
        print '\t',np.abs(_y0),N_px 
        print '\t',_amplitude0, 0., 100.*amp1
        print '\t',_amplitude1, 0., 100.*amp2
        print '\t',_amplitude2, 0., 100.*amp3
        print '\t',_amplitude3, 0., 100.*amp4
        print '\t',_r0_lens, 1.35/float(FM.hdulist[0].header['PIXSCALE']), 1.65/float(FM.hdulist[0].header['PIXSCALE']) 
        print '\t',_q_lens, 1.002, 1.006 
        print '\t',_posang_lens, 171.4, 178.2
        print '\t',_alpha, 0., 100.*alpha
        print '\t',_beta,  0., 100.*beta 
        print '\t',_N_bkg, 0., 100.*N_bkg_guess
        print 'lnp is -infinity'
        '''
        return -np.inf

  def lnprob(theta, obj, N_pix, flip):
    lp = lnprior(theta)
    #print lp, theta
    if not np.isfinite(lp):
        return -np.inf
    #mll = obj.moffat_chi_sq(theta, N_pix, flip, x_images, y_images)
    mll = obj.moffat_chi_sq_w_lens(theta, N_pix, flip, x_images, y_images, x_lens, y_lens)
    #print mll
    if not np.isfinite(mll):
        return -np.inf
    #print '%1.2e\t%1.2e'%(lp,mll), theta
    return lp - mll



  #################################################
  #################################################
  #################################################

  # RUN MC!
  x0,y0,amp1,amp2,amp3,amp4, amp_lens, alpha, beta, N_bkg = FM.qsr_min_parms
  theta = [x0,y0,amp1,amp2,amp3,amp4, amp_lens, r0_lens, q_lens, posang_lens, alpha, beta, N_bkg]
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
  labels = [ "x0", "y0", "amp1", "amp2", "amp3", "amp4", "amp_lens", "r0_lens", "q_lens", "posang_lens", "alpha", "beta", "N_bkg"]
  fig= corner.corner(samples, labels=labels)
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
  for k in range(0,n_iterations/1000):
      sampler.run_mcmc(pos, 1000)
      print ''
      print k+1, (sampler.chain).shape, datetime.datetime.now()
      samples = sampler.chain[:, int(0.5*float((k+1)*100)):, :].reshape((-1, ndim))
      print("\tMean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
      parms = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))

      #print parms.shape
      #print parms[k][0]

      for k in range(0,len(parms)):
        print '%s\t%+1.2e\t%+1.2e\t%+1.2e\t%+1.2e'%(labels[k], theta[k], parms[k][0],parms[k][1],parms[k][2])

      fig= corner.corner(samples, labels=labels)
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

  '''
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
  '''
def flux2magnitude(flux):
    #for k in range(0,len(flux)):
    #    if(flux[k]<=0.):
    #        print 'bad',k,flux[k]
    return -2.5*np.log10(flux)

def magnitude2flux(mag):
    return np.power(10.,-mag/2.5)

def magErr2fluxErr(magVal, magErr):
    flux = magnitude2flux(magVal)
    return flux*(1.-np.power(10.,-magErr/2.5))

def fluxErr2magErr(fluxVal, fluxErr):
    return -2.5*np.log10(1-fluxErr/fluxVal)
    




# UNUSED FUNCTIONS #
################################################################
################################################################
################################################################

################################################################
################################################################
################################################################

################################################################
################################################################
################################################################

################################################################
################################################################
################################################################

################################################################
################################################################
################################################################

################################################################

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

################################################################

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

################################################################

def triple_gaussian((x, y), amp1, sigx1, sigy1, ang1, xo1, yo1, amp2, sigx2, sigy2, ang2, xo2, yo2, amp3, sigx3, sigy3, ang3, xo3, yo3):
	g1 = gaussian((x, y), amp1, sigx1, sigy1, ang1, xo1, yo1)
	g2 = gaussian((x, y), amp2, sigx2, sigy2, ang2, xo2, yo2)
	g3 = gaussian((x, y), amp3, sigx3, sigy3, ang3, xo3, yo3)
	return g1 + g2 + g3

################################################################

def triple_gauss_chi2(parms, image, sig):
	xg, yg = np.mgrid[:len(image), :len(image)]
	model = triple_gaussian( (xg,yg,), *parms)
	chi = (image-model)/sig
	print np.sum(chi*chi), parms
	return np.sum(chi*chi)

################################################################

def fitter(image, p0):
    sig = np.sqrt(image)
    results = minimize(triple_gauss_chi2, p0, args=(image, sig), method='Nelder-Mead')
    print 'finished running minimize'
    # return the best fit parameters
    print results
    return results.x
#get_background_counts()
#exit()
################################################################



def twoD_Moffat((x, y), amplitude, alpha, beta, xo, yo, offset, pixel_integration = True):
    #print 'len(x), len(y)', len(x), len(y)
    if(pixel_integration==True):
        xo = float(xo - (len(x)- 1) / 2) 
        yo = float(yo - (len(y)- 1) / 2)
        PSF_model = psf.moffat_kernel(xo, yo, alpha, beta, nx=len(x), ny=len(y))
        PSF_model /= np.max(PSF_model)
        PSF_model = np.rot90(PSF_model)
        PSF_model = np.flipud(PSF_model)    
        m = amplitude*PSF_model + offset
        if(alpha<0.): m+=1.e9
        if(beta<0.): m+=1.e9
        return m.ravel() 
    if(pixel_integration==False):
        xo = float(xo)
        yo = float(yo)    
        a = (beta-1.)/(np.pi*alpha**2)
        m = offset + amplitude*( 1. + ((x-xo)**2 + (y-yo)**2) / (alpha**2))**(-beta)
        if(alpha<0.): m+=1.e9
        if(beta<0.): m+=1.e9
        return m.ravel()

################################################################

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

################################################################

def twoD_Moffat_proj(x, amplitude, alpha, beta, xo, yo, offset):
  a = (beta-1.)/(np.pi*alpha**2)
  m = offset + amplitude * ( 1. + (x**2) / (2.*alpha**2))**(-beta)
  #print np.array(g)
  #print g
  return m.ravel()
  #return np.array(g)    

################################################################

def moffat_analytical_integral(amplitude, alpha, beta):
  return amplitude * 2. * np.pi * alpha**2 / (beta-1.)

################################################################

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


################################################################

def twoD_Gaussian_simple((x, y), amplitude, xo, yo, sigma, offset):
  xo = float(xo)
  yo = float(yo)    
  a = 1./(2*sigma**2)
  g = offset + amplitude*np.exp( - ((x-xo)**2 + (y-yo)**2) / (2.*sigma**2))
  #print np.array(g)
  #print g
  return g.ravel()
  #return np.array(g)  
  
################################################################

def fit_gaussian(image, verbose=False):

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

################################################################

def radial_profile(image, x0,y0):
  r=[]
  amp=[]
  for i in range(0,len(image)):
      for j in range(0,len(image)):
	r.append(np.sqrt((float(i)-x0)**2 + (float(j)-y0)**2))
	amp.append(image[i][j])
  return np.array(r), np.array(amp)

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

################################################################

def moffat_fwhm(alpha,beta):
	return alpha*2.*np.sqrt(2.**(1/beta)-1.)

################################################################

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


################################################################

# The functions below were part of the SourceImage class
def twoD_Moffat_chi(self, theta, readnoise):
    model = twoD_Moffat((self.xg, self.yg), *theta).reshape(len(self.image),len(self.image))
    return (self.image-model)/np.sqrt(self.image + readnoise**2)

################################################################

def twoD_Moffat_chi_sq(self, theta ):
    return np.sum(self.twoD_Moffat_chi(theta, self.FM.readnoise)**2)

################################################################

def fit_moffat(self, readnoise, verbose=False):
    if(verbose): print 'fitting'
    self.x_max, self.y_max = np.unravel_index(self.image.argmax(), self.image.shape)
    # estimate fwhm
    pk = np.max(self.image)
    fwhm=1.
    print 'pk, fwhm', pk, fwhm
    for dx in range(0,15):
        val_a = self.image[self.x_max + dx, self.y_max]
        val_b = self.image[self.x_max + dx+1, self.y_max]
        if(val_a>=pk/2. and val_b<=pk/2.):
	     fwhm = dx+0.5
	     break
    print 'crude fwhm', fwhm
    print 'x_max, y_max', self.x_max, self.y_max
    x0_guess = self.x_max
    y0_guess = self.y_max
    print 'x_max, y_max', self.x_max, self.y_max
    print 'x0_guess, y0_guess', x0_guess, y0_guess
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
    guess_beta = 2.
    guess_alpha = fwhm/2./(2.0 ** (1.0 / guess_beta) - 1.0) ** 0.5  
    print 'guess_alpha', guess_alpha

    self.initial_guess_moffat = [self.image[self.x_max,self.y_max]+background_count_guess, guess_alpha, guess_beta, x0_guess, y0_guess, background_count_guess]
    print self.initial_guess_moffat
    if(verbose): print 'initial guess', self.initial_guess_moffat
    if(verbose): print len(self.image)
    x=np.arange(0.,len(self.image))
    y=np.arange(0.,len(self.image))
    self.xg, self.yg = np.mgrid[:len(self.image), :len(self.image)]
    #popt, pcov = opt.curve_fit(twoD_Gaussian, (xg, yg), image.ravel(), p0=initial_guess)
    #popt, pcov = opt.curve_fit(twoD_Gaussian_simple, (xg, yg), image.ravel(), p0=initial_guess_simple)
    #try:
    #results = opt.minimize(self.twoD_Moffat_chi_sq, self.initial_guess_moffat, method='Nelder-Mead')
    #popt = results.x
    #exit()	
    popt, pcov = opt.curve_fit(twoD_Moffat, (self.xg, self.yg), self.image.ravel(), p0=self.initial_guess_moffat)
    print 'popt',popt
    print 'pcov',np.sqrt(pcov[0][0]), np.sqrt(pcov[0][0])/popt[0], np.linalg.det(pcov)
    for i in range(0,len(popt)):
       for j in range(0,len(popt)):
           pcov[i][j]/=popt[i]*popt[j]
    print 'stat_cut?', np.power(np.sqrt(np.linalg.det(pcov)), 1./len(popt))
    print 'stat_cut?', np.sqrt(np.trace(pcov)/len(popt))
    if ( np.sqrt(np.trace(pcov)/len(popt)) > 0.15):
        self.fit_ok = False
    else:
        self.fit_ok = True
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
    #print 'shapes', self.image.shape, self.moffat_fit_image.shape
    print 'readnoise', readnoise
    self.moffat_chi = (self.image-self.moffat_fit_image)/np.sqrt(self.image + readnoise**2)
    print '\tmoffat estimate chi_sq', np.sum(self.moffat_chi**2)/len(self.moffat_chi)
    return popt

################################################################

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

################################################################

def plot_moffat_residual_2D(self):
    plt.imshow(self.image-self.moffat_fit_image, cmap='gray', interpolation='none')
    plt.title('ra: %1.4f dec: %1.4f'%(self.ra,self.dec))
    plt.colorbar()

################################################################

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

################################################################

def plot_radial_profile_residuals(self):
    plt.plot(self.rad,self.moffat_chi,'b.')
    plt.title('ra: %1.4f dec: %1.4f'%(self.ra,self.dec))
    plt.ylim(-15.,15.)
    plt.xlabel('Pixel Distance from Fitted Peak')
    plt.ylabel(r'$\chi_p$', fontsize=20)

################################################################

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

################################################################

def moffat_chi_vals(self,theta,N_pix, flip, x_images, y_images):
    x0,y0,amp0,amp1,amp2,amp3, alpha, beta, N_bkg = theta
    if(amp0<0 or amp1<0 or amp2<0 or amp3<0):
	    return np.inf

    model = self.quad_image_model(x0,y0,x_images,y_images, amp0,amp1,amp2,amp3, alpha, beta, N_bkg, N_pix, flip)
    chi = (self.image - model)/np.sqrt(self.image+self.FM.readnoise**2)
    #print np.sum(chi*chi)
    #print x0,y0, N_pix/2
    # NO NEGATIVE AMPLITUDES ALLOWED!
    #print 'moffat_chi_vals: %1.5e %1.2e %1.2e %1.2e %1.2e %1.2e'%(np.sum(chi**2), amp0, amp1, amp2, amp3, N_bkg)
    # THE IMAGE CORRECTION HAS TO BE WITHIN THE BOUNDS OF THE IMAGE!
    if(x0>N_pix/2 or x0<-N_pix/2 or y0>N_pix/2 or y0<-N_pix/2):
	    return np.inf
    return chi

################################################################

def moffat_chi_sq(self, theta, N_pix, flip, x_images, y_images):
	chisq = np.sum((self.moffat_chi_vals(theta, N_pix, flip, x_images, y_images))**2)
	#print '%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e'%(chisq, theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], theta[6], theta[7], theta[8]) 
	return chisq

################################################################

def moffat_chi_vals_w_lens(self,theta,N_pix, flip, x_images, y_images, x_lens, y_lens):
    x0,y0,amp0,amp1,amp2,amp3, amp_lens, r0_lens, q_lens, posang_lens, alpha, beta, N_bkg = theta
    if(amp0<0 or amp1<0 or amp2<0 or amp3<0 or amp_lens<0. or r0_lens<=0. or q_lens<0. or posang_lens<0.):
	    return np.inf

    model = self.quad_image_model(x0,y0,x_images,y_images, amp0,amp1,amp2,amp3, alpha, beta, N_bkg, N_pix, flip, x_lens, y_lens, amp_lens, r0_lens, q_lens, posang_lens)
    chi = (self.image - model)/np.sqrt(self.image+self.FM.readnoise**2)
    #print np.sum(chi*chi)
    #print x0,y0, N_pix/2
    # NO NEGATIVE AMPLITUDES ALLOWED!
    #print 'deVaucouleurs Model Peak', np.max(m)

    m = deVaucouleurs_model(x0, y0, amp_lens, r0_lens, alpha, beta, q = q_lens, posang=posang_lens, npx=31)
    #print 'moffat_chi_vals_w_lens: %1.5e %1.2e %1.2e %1.2e %1.2e %1.2e %1.2e %1.2e %1.1f %1.1f'%(np.sum(chi**2), amp0, amp1, amp2, amp3, N_bkg, amp_lens, r0_lens, posang_lens, np.max(m))
    # THE IMAGE CORRECTION HAS TO BE WITHIN THE BOUNDS OF THE IMAGE!
    if(x0>N_pix/2 or x0<-N_pix/2 or y0>N_pix/2 or y0<-N_pix/2):
	    return np.inf
    return chi

################################################################

def moffat_chi_sq_w_lens(self, theta, N_pix, flip, x_images, y_images, x_lens, y_lens):
	chisq = np.sum((self.moffat_chi_vals_w_lens(theta, N_pix, flip, x_images, y_images, x_lens, y_lens))**2)
	#print '%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e'%(chisq, theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], theta[6], theta[7], theta[8]) 
	return chisq

################################################################

def moffat_chi_vals_w_lens_amp_only(self,theta,N_pix, flip, x_images, y_images, x_lens, y_lens, r0_lens, q_lens, posang_lens):
    x0,y0,amp0,amp1,amp2,amp3, amp_lens, alpha, beta, N_bkg = theta
    if(amp0<0 or amp1<0 or amp2<0 or amp3<0 or amp_lens<0. or r0_lens<=0. or q_lens<0. or posang_lens<0.):
	    return np.inf

    model = self.quad_image_model(x0,y0,x_images,y_images, amp0,amp1,amp2,amp3, alpha, beta, N_bkg, N_pix, flip, x_lens, y_lens, amp_lens, r0_lens, q_lens, posang_lens)
    chi = (self.image - model)/np.sqrt(self.image+self.FM.readnoise**2)
    #print np.sum(chi*chi)
    #print x0,y0, N_pix/2
    # NO NEGATIVE AMPLITUDES ALLOWED!
    #print 'deVaucouleurs Model Peak', np.max(m)

    m = deVaucouleurs_model(x0, y0, amp_lens, r0_lens, alpha, beta, q = q_lens, posang=posang_lens, npx=31)
    #print 'moffat_chi_vals_w_lens: %1.5e %1.2e %1.2e %1.2e %1.2e %1.2e %1.2e %1.2e %1.1f %1.1f'%(np.sum(chi**2), amp0, amp1, amp2, amp3, N_bkg, amp_lens, r0_lens, posang_lens, np.max(m))
    # THE IMAGE CORRECTION HAS TO BE WITHIN THE BOUNDS OF THE IMAGE!
    if(x0>N_pix/2 or x0<-N_pix/2 or y0>N_pix/2 or y0<-N_pix/2):
	    return np.inf
    return chi

################################################################

def moffat_chi_sq_w_lens_amp_only(self, theta, N_pix, flip, x_images, y_images, x_lens, y_lens, r0_lens, q_lens, posang_lens):
	chisq = np.sum((self.moffat_chi_vals_w_lens_amp_only(theta, N_pix, flip, x_images, y_images, x_lens, y_lens, r0_lens, q_lens, posang_lens))**2)
	#print '%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e'%(chisq, theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], theta[6], theta[7], theta[8]) 
	return chisq

'''
def moffat_chi_vals2(self,theta,alpha, beta,N_pix, flip):
    x0,y0,amp0,amp1,amp2,amp3, N_bkg = theta
    model = self.quad_image_model(x0,y0,amp0,amp1,amp2,amp3, alpha, beta, N_bkg, N_pix, flip)

    chi = (self.image - model)/np.sqrt(self.image+self.FM.readnoise**2)
    #print 'self.FM.readnoise', self.FM.readnoise
    #print np.sum(chi*chi)
    return chi
def moffat_chi_sq2(self, theta, alpha, beta, N_pix, flip):
	return np.sum((self.moffat_chi_vals2(theta, alpha, beta, N_pix, flip))**2)
'''


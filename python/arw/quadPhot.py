#!/usr/bin/env python

'''
2015 March 15
Andrew Romero-Wolf
Estimation of quadruply lensed quasar images with LCOGT images.
'''

import matplotlib
matplotlib.use('Agg') 
from AnalysisLibrary import * 
import sys
import argparse
matplotlib.rcParams['figure.facecolor']='white'

if __name__ == "__main__":

	parser=argparse.ArgumentParser(description='quadPhot routine to calculate crowded field quadruply lensed quasar photometry')
	parser.add_argument("-d","--dataDir", default='/nisushome/data2/romerowo/lcogt_data/he045-1223_wcs_corrected', help="directory with FITS file images",type=str)
	parser.add_argument("-i","--inputFile", default = 'lsc1m009-fl03-20141217-0042-e90.fits', help="directory with FITS file images",type=str)
	parser.add_argument("-o","--outputFileTag", default = 'out', help="directory with FITS file images",type=str)
	parser.add_argument("-p","--plots", default = 1, help="level of plot production",type=int)
	parser.add_argument("-e","--emcee_level", default = 0, help="level of emcee fitting (0 none, 1 quad image, 2 quad image w/ photometry, 3 full set images w/ photometry",type=int)

	#KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
	for i, arg in enumerate(sys.argv):
	  if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

	args=parser.parse_args()
	inputFile = args.dataDir+'/'+args.inputFile
	print '\n############################################################'
	print 'quadPhot inputFile', inputFile.split('/')[-1]
	print 'quadPhot outputFile',args.outputFileTag
	print '############################################################\n'
	print args
	 
	# DECIMAL RA and DEC VALUES OF HE0435-1223
	ra_qsr = (4.+38./60.+14.9/60./60.)/24*360 
	dec_qsr = -12. - (17./60 +14.4/60./60.)

	npxls = 31
	# INITIALIZE THE CUSTOMIZED FITS FILE MANAGER
	FM = FITSmanager(inputFile)

	# PLOT QSR IMAGE
	FM.plot_image(ra_qsr, dec_qsr, Npx=npxls, out = args.outputFileTag+'_qsr')

	# ESTIMATE THE BACKGROUND AND READ NOISE
	FM.estimate_read_noise(display=args.plots, out = args.outputFileTag+'_readnoise')

	# APASS PHOTOMETRY
	APASS_table = ascii.read('../../data/HE0435_LCOGT/APASS_0438_list.csv')
	APASS_rejects = [9, 21, 22] # 9 and 21 seem to be at the edge of the field of view. 22 seems close to the edge, sometimes we don't catch it.
	ZP_mean, ZP_wrms, ZP_rms, alpha_mean, beta_mean, alpha_wrms, beta_wrms, alpha_beta_corr = APASS_zero_points(FM, APASS_table, APASS_rejects, FM.readnoise, display=args.plots, out = args.outputFileTag+'_APASS')
	print '\tZero Points',ZP_mean, ZP_wrms
	print '\tSeeing alpha, beta, err_alpha, err_beta, corr_alpha_beta',alpha_mean, beta_mean, alpha_wrms, beta_wrms, alpha_beta_corr

	# FROM CASTLES
	ra_images=np.array([0.,-1.476,-2.467,-0.939])
	dec_images=np.array([0.,0.553, -0.603, -1.614])
	ra_images -= np.mean(ra_images)
	dec_images -= np.mean(dec_images)
	# FROM Wisotzki02
	ra_imagesW=[0.,-1.483,-2.488,-0.951]
	dec_imagesW=[0.,0.567, -0.589, -1.620]
	ra_imagesW -= np.mean(ra_imagesW)
	dec_imagesW -= np.mean(dec_imagesW)
	print ra_images - ra_imagesW
	print dec_images - dec_imagesW

	print 'PIXSCALE', FM.hdulist[0].header['PIXSCALE']
	print 'ra images in arcsec', ra_images
	print 'dec images in arcsec', dec_images
	m1, me1, m2, me2, m3, me3, m4, me4, chiSq, maxChi = quadFit(FM, ra_qsr, dec_qsr, ra_images, dec_images, ZP_mean, ZP_wrms, alpha_mean, beta_mean, npxls, outputFileTag=args.outputFileTag)

	npz_out = args.outputFileTag + '_results.npz'
	readnoise = FM.readnoise
	APASS_alpha = alpha_mean
	APASS_beta  = beta_mean
	APASS_alpha_err = alpha_wrms
	APASS_beta_err  = beta_wrms
	APASS_alpha_beta_corr  = alpha_beta_corr

	outFileTag = args.outputFileTag
	mjd_obs = float(FM.hdulist[0].header['MJD-OBS'])
	# INCLUDE ALL APASS FIT RESULTS, alpha, beta, chiSq, maxChi
	# SAVE RESULTS TO AN NPZ FILE
	filter = FM.hdulist[0].header['FILTER']
	np.savez(npz_out, 
	inputFile = inputFile,
	outFileTag = outFileTag, 
	mjd_obs = mjd_obs,
	readnoise = readnoise, 
	ZP_mean = ZP_mean, 
	ZP_wrms = ZP_wrms, 
	ZP_rms = ZP_rms, 
	APASS_alpha = APASS_alpha,
	APASS_beta = APASS_beta,
	APASS_alpha_err = APASS_alpha_err,
	APASS_beta_err = APASS_beta_err,
	APASS_alpha_beta_corr = APASS_alpha_beta_corr,
	filter = filter,
	m1 = m1,
	me1 = me1,
	m2 = m2,
	me2 = me2,
	m3 = m3,
	me3 = me3,
	m4 = m4,
	me4 = me4,
	chiSq = chiSq,
	maxChi = maxChi
	)

	if(args.emcee_level==1):
		emceeQuadFit(FM, ra_qsr, dec_qsr, ZP_mean, ZP_wrms, alpha_mean, beta_mean, alpha_wrms, beta_wrms, alpha_beta_corr, m1, m2, m3, m4, npxls, ra_images, dec_images,  outputFileTag=args.outputFileTag)
	if(args.emcee_level>1):
		print 'YOU HAVE SELECTED A LEVEL OF EMCEE FITTING OF %d'%(args.emcee_level)
		print 'THIS CAPABILITY HAS NOT YET BEEN DEVELOPED'
		print 'THANK YOU FOR YOUR PATIENCE'
	plt.show()



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
    parser.add_argument("-d","--dataDir", default='/halo_nobackup/lenssim/romerowo/lcogt_data/Dec_2015_100_hours', help="directory with FITS file images",type=str)
    parser.add_argument("-i","--inputFile", default = 'coj1m003-kb71-20151214-0083-e91.fits', help="directory with FITS file images",type=str)
    parser.add_argument("-nmax","--nmax", default = 12, help="Maximum n mode for shapelet fitting",type=int)
    parser.add_argument("-mmax","--mmax", default = 3, help="Maximum m mode for shapelet fitting",type=int)
    parser.add_argument("-o","--outputFileTag", default = 'out', help="directory with FITS file images",type=str)
    parser.add_argument("-p","--plots", default = 1, help="level of plot production",type=int)
    parser.add_argument("-sl","--star_list", default = 'proofStars.txt', help="list of stars to test the PSF fit on",type=str)
    parser.add_argument("-e","--emcee_level", default = 0, help="level of emcee fitting (0 none, 1 quad image, 2 quad image w/ photometry, 3 full set images w/ photometry",type=int)

    #KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
    for i, arg in enumerate(sys.argv):
      if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

    args=parser.parse_args()
    inputFile = args.dataDir+'/'+args.inputFile
    print '\n############################################################'
    print 'quadPhot inputFile ', inputFile.split('/')[-1]
    print 'quadPhot outputFile',args.outputFileTag
    print '############################################################\n'
    # print args
    # print 
    print '\nInputs Parameters and Values'
    for arg in vars(args):
        print '\t', arg.ljust(14), '= ', getattr(args, arg)
    # DECIMAL RA and DEC VALUES OF HE0435-1223
    ra_qsr = (4.+38./60.+14.9/60./60.)/24*360 
    dec_qsr = -12. - (17./60 +14.4/60./60.)

    npxls = 31
    # INITIALIZE THE CUSTOMIZED FITS FILE MANAGER
    FM = FITSmanager(inputFile, args.nmax, args.mmax)

    # PLOT QSR IMAGE
    # FM.plot_image(ra_qsr, dec_qsr, Npx=npxls, out = args.outputFileTag+'_qsr')

    # ESTIMATE THE BACKGROUND AND READ NOISE
    FM.estimate_read_noise(display=args.plots, out = args.outputFileTag+'_readnoise')

    # APASS PHOTOMETRY
    APASS_table = ascii.read('/halo_nobackup/lenssim/romerowo/lcolens/python/arw/APASS_0438_list.csv')
    #APASS_rejects = [9, 21, 22] # 9 and 21 seem to be at the edge of the field of view. 22 seems close to the edge, sometimes we don't catch it.
    APASS_rejects = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17, 19, 21, 22] # 1-17 are too bright. 18 and 19 are repeats, 9 and 21 seem to be at the edge of the field of view. 22 seems close to the edge, sometimes we don't catch it.
    APASS_index_list, beta_mean, beta_unc, shapelet_coeffs, APASS_S_CCD_list, APASS_sig_S_CCD_list, APASS_chi_sq, APASS_max_chi = APASS_zero_points(FM, APASS_table, APASS_rejects, FM.readnoise, display=args.plots, out = args.outputFileTag+'_APASS')


    # FROM CASTLES
    ra_images=np.array([0.,-1.476,-2.467,-0.939])
    dec_images=np.array([0.,0.553, -0.603, -1.614])
    ra_lensgal  = -1.165
    dec_lensgal = -0.573
    ra_offset  = np.mean(ra_images)
    dec_offset = np.mean(dec_images)
    ra_images -= ra_offset
    dec_images -= dec_offset
    ra_lensgal  -= ra_offset
    dec_lensgal -= dec_offset
    '''
    # FROM Wisotzki02
    ra_imagesW=[0.,-1.483,-2.488,-0.951]
    dec_imagesW=[0.,0.567, -0.589, -1.620]
    ra_imagesW -= np.mean(ra_imagesW)
    dec_imagesW -= np.mean(dec_imagesW)
    ra_lensgalW  = -1.165
    dec_lensgalW = -0.573
    ra_offsetW  = np.mean(ra_imagesW)
    dec_offsetW = np.mean(dec_imagesW)
    ra_imagesW -= ra_offsetW
    dec_imagesW -= dec_offsetW
    ra_lensgalW  -= ra_offsetW
    dec_lensgalW -= dec_offsetW
    #print ra_images - ra_imagesW
    #print dec_images - dec_imagesW
    '''

    #popt_ng, pcov_ng, chisq_ng, max_chi_ng = quadFit(FM, ra_qsr, dec_qsr, ra_images, dec_images, ra_lensgal, dec_lensgal, beta_mean, shapelet_coeffs, npxls, galFit=False, display=args.plots,  outputFileTag=args.outputFileTag, emcee_level = args.emcee_level)


    popt_wg, pcov_wg, chisq_wg, max_chi_wg = quadFit(FM, ra_qsr, dec_qsr, ra_images, dec_images, ra_lensgal, dec_lensgal, beta_mean, shapelet_coeffs, npxls, galFit=True, display=args.plots,  outputFileTag=args.outputFileTag, emcee_level = args.emcee_level)

    star_ra, star_dec = readStarList(args.star_list)
    star_index_list, star_chi_sq, star_max_chi, star_S_CCD, star_S_CCD_unc = starFit(FM, star_ra, star_dec, beta_mean, shapelet_coeffs, N_px=npxls, display = args.plots, outputFileTag=args.outputFileTag)

    #exit()

    filter   = FM.hdulist[0].header['FILTER']
    mjd_obs  = float(FM.hdulist[0].header['MJD-OBS'])
    npz_out  = args.outputFileTag + '_results.npz'

    np.savez(npz_out,
	inputFile       = inputFile,
	outFileTag      = args.outputFileTag,
	mjd_obs         = mjd_obs,
	readnoise       = FM.readnoise,
	beta_mean       = beta_mean,
	beta_unc        = beta_unc,
	shapelet_coeffs = shapelet_coeffs,
	nmax            = args.nmax,
	mmax            = args.mmax,
	filter          = filter,
	APASS_index_list= APASS_index_list,
	APASS_S_CCD     = APASS_S_CCD_list,
	APASS_sig_S_CCD = APASS_sig_S_CCD_list,
	APASS_chi_sq    = APASS_chi_sq, 
	APASS_max_chi   = APASS_max_chi,
    star_index_list = star_index_list,
	star_S_CCD      = star_S_CCD,
	star_S_CCD_unc  = star_S_CCD_unc,
	star_chi_sq     = star_chi_sq,
	star_max_chi    = star_max_chi,
    qsr_wg_parms    = popt_wg, 
    qsr_wg_covar    = pcov_wg,
    qsr_chisq       = chisq_wg, 
    qsr_max_chi     = max_chi_wg
    )
    print 'quadPhot Complete'
    #'''
    ######
    exit()
    ######

    if(args.emcee_level==1):
	    emceeQuadFit(FM, ra_qsr, dec_qsr, ZP_mean, ZP_wrms, alpha_mean, beta_mean, alpha_wrms, beta_wrms, alpha_beta_corr, m1, m2, m3, m4, npxls, ra_images, dec_images,  outputFileTag=args.outputFileTag)
    if(args.emcee_level==2):
	    emceeQuadFitLensGal(FM, ra_qsr, dec_qsr, ZP_mean, ZP_wrms, alpha_mean, beta_mean, alpha_wrms, beta_wrms, alpha_beta_corr, m1, m2, m3, m4, npxls, ra_images, dec_images, ra_lensgal, dec_lensgal, outputFileTag=args.outputFileTag)
    if(args.emcee_level>2):
	    print 'YOU HAVE SELECTED A LEVEL OF EMCEE FITTING OF %d'%(args.emcee_level)
	    print 'THIS CAPABILITY HAS NOT YET BEEN DEVELOPED'
	    print 'THANK YOU FOR YOUR PATIENCE'
    plt.show()



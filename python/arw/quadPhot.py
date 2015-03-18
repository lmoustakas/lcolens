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

if __name__ == "__main__":
    
    parser=argparse.ArgumentParser(description='quadPhot routine to calculate crowded field quadruply lensed quasar photometry')
    parser.add_argument("-d","--dataDir", default='/data2/romerowo/lcogt_data/he045-1223_wcs_corrected', help="directory with FITS file images",type=str)
    parser.add_argument("-i","--inputFile", default = 'lsc1m009-fl03-20141217-0042-e90.fits', help="directory with FITS file images",type=str)
    parser.add_argument("-o","--outputFileTag", default = 'out', help="directory with FITS file images",type=str)
    parser.add_argument("-p","--plots", default = 1, help="level of plot production",type=int)
    
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
    ZP_mean, ZP_rms, alpha_mean, beta_mean = APASS_zero_points(FM, APASS_table, APASS_rejects, FM.readnoise, display=args.plots, out = args.outputFileTag+'_APASS')
    print ZP_mean, ZP_rms, alpha_mean, beta_mean

    m1, me1, m2, me2, m3, me3, m4, me4, chiSq, maxChi = quadFit(FM, ra_qsr, dec_qsr, ZP_mean, ZP_rms, alpha_mean, beta_mean, npxls, outputFileTag=args.outputFileTag)
    
    npz_out = args.outputFileTag + '_results.npz'
    readnoise = FM.readnoise
    APASS_alpha = alpha_mean
    APASS_beta  = beta_mean

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
	ZP_rms = ZP_rms, 
	APASS_alpha = APASS_alpha,
	APASS_beta = APASS_beta,
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
    # plt.show()



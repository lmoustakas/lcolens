#!/usr/bin/env python

'''
2015 April 29
Andrew Romero-Wolf
Estimation of quadruply lensed quasar images with LCOGT images.
'''
# 2015, 04/29: file created as a cleaned up version of lcolens/python/arw/quadPhot.py

import matplotlib
matplotlib.use('Agg') 
import astropy

#import analysis_library as AL 
#import sys
#import argparse

if __name__ == "__main__":
    
    parser=argparse.ArgumentParser(description='quadPhot routine to calculate crowded field quadruply lensed quasar photometry')
    parser.add_argument("-d","--dataDir", default='/data2/romerowo/lcogt_data/he045-1223_wcs_corrected', help="directory with FITS file images",type=str)
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
 
    # Use as much astropy as possible.
	# hoping that PSF development used for this project can be pushed back into astropy.
	# 

    # 1. Define global parameters to be used in the analysis (use a class for this).
    #	  i. define frame size in arcseconds that will be used for source images 

    # 2. Read FITS file
    #	   i. read parameters that will be used throughout the analysis (arcseconds per pixel, etc.)
    #	  ii. estimate read noise as a separate routine.

    # 3. Derive Zero points
    #      i. define class that loads APASS star parameters
    #     ii. apply routines fitstars.py to fit profiles to reference stars in FITS file
    #    iii. add airmass terms. include extinction coefficient to tanslate between the magnitude above the atmosphere rather than below the atmosphere.
    #     iv. color terms.  (see Sloan Photometry paper).
    #	   v. determine total light from fitted profile. 
    #     vi. convert to flux and magnitudes
    #    vii. Write a note on what is being done.
    # 	viii. derive zero points (what fitting approach will we use here?) wMean, wRMS.

    # 4. Quadruple Image Fitting
    # 	  i. define template for the 4 images. Template registration (HST to LCOGT). Either use measured positions or use the one from HST. The HST registration need to be done on an image by image basis. Do this once and derive the uncertainties. Use it as a prior for the photometry.
    #    ii. apply an adaptation of routines in fitstars.py to the quad image fitter.



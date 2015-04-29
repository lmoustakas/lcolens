#!/usr/bin/env python

'''
2015 April 29
Andrew Romero-Wolf
Estimation of quadruply lensed quasar images with LCOGT images.
'''
# 2015, 04/29: file created as a cleaned up version of lcolens/python/arw/quadPhot.py

import matplotlib
matplotlib.use('Agg') 
import analysis_library as AL 
import sys
import argparse

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
 
    # 1. Define global parameters to be used in the analysis.
    #	  i. define frame size in arcseconds that will be used for source images 

    # 2. Read FITS file
    #	   i. read parameters that will be used throughout the analysis (arcseconds per pixel, etc.)
    #	  ii. estimate read noise

    # 3. Derive Zero points
    #      i. define class that loads APASS star parameters
    #     ii. apply routines fitstars.py to fit profiles to reference stars in FITS file
    #	 iii. determine total light from fitted profile. 
    #     iv. convert to flux and magnitudes
    # 	   v. derive zero points (what fitting approach will we use here?)

    # 4. Quadruple Image Fitting
    # 	  i. define template for the 4 images
    #    ii. apply an adaptation of routines in fitstars.py to the quad image fitter.



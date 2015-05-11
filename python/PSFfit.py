#!/usr/bin/env python

'''
2015 May 06
Andrew Romero-Wolf
PSF fitting to extract parameters and their uncertainties
'''

import matplotlib
#matplotlib.use('Agg')
import pylab 
import sys
import astropy.io.ascii
import astropy.io.fits
import astropy.wcs
import argparse
import analysis_library as AL
import numpy as np

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
    print 'PSFfit inputFile', inputFile.split('/')[-1]
    print 'PSFfit outputFile',args.outputFileTag
    print '############################################################\n'
    print args
 
    #############
    # Code flow #
    #############

    # 1. Set up global variables
    # 2. Read APASS source catalogue
    # 3. Read FITS image file
    # 4. Estimate image noise parameters
    # 5. PSF fitting
    #    5.a. First step, read single APASS stars and fits them individually. Extract the weighted mean and RMS of their parameters.
    #    5.b. Second step, fit all stars simultaneously using a minimizer.
    #    5.c. Third step, fit all stars simultaneously using emcee. Return a fit to the distribution of PSF parameters. 
    
    # 1. Set up global variables
    global_parameters = AL.PSFglobals(args.dataDir, args.inputFile, args.outputFileTag, args.plots) 

    # 2. Read APASS source catalogue
    APASS_table = astropy.io.ascii.read('../data/HE0435_LCOGT/APASS_0438_list.csv')
    print APASS_table

    # 3. Read FITS image file
    hdulist = astropy.io.fits.open(inputFile)
    FITSdata = hdulist[0]
    print hdulist.info()
    '''
    count=0
    for k in FITSdata.header.keys():
	count+=1
        spaces = ' '*12
	print '%s%s%s'%(k,spaces[len(k):], FITSdata.header[k])
	if(count>=100): break
    print FITSdata.header['MJD-OBS']
    '''

    # 4. Estimate image noise parameters
    # second level astrometric alignment pixel registration
    img_arcsec_width=16.
    Npx = np.ceil(img_arcsec_width/FITSdata.header['PIXSCALE'])
    if(Npx%2==0): Npx+=1
    #for k in range(0,3):
    for k in range(0,len(APASS_table)):
        print k, APASS_table['radeg'][k], APASS_table['decdeg'][k]
        AL.plot_subimage(FITSdata, APASS_table['radeg'][k], APASS_table['decdeg'][k], Npx)
    pylab.show()
    print Npx
    # px_ra, px_dec, counts = getSubImage(x_px, y_px, nx, ny)
   





    


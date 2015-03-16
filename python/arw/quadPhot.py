#!/usr/bin/env python

'''
2015 March 15
Andrew Romero-Wolf
Estimation of quadruply lensed quasar images with LCOGT images.
'''

from AnalysisLibrary import * 
import sys
import argparse
import matplotlib

if __name__ == "__main__":
    
    parser=argparse.ArgumentParser(description='quadPhot routine to calculate crowded field quadruply lensed quasar photometry')
    parser.add_argument("-d","--dataDir", default='/data2/romerowo/lcogt_data/he045-1223_wcs_corrected', help="directory with FITS file images",type=str)
    parser.add_argument("-i","--inputFile", default = 'lsc1m009-fl03-20141217-0042-e90.fits', help="directory with FITS file images",type=str)
    parser.add_argument("-o","--outputFileTag", default = 'out', help="directory with FITS file images",type=str)
    
    #KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
    for i, arg in enumerate(sys.argv):
      if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

    args=parser.parse_args()
    print args
    inputFile = args.dataDir+'/'+args.inputFile
    print 'inputFile', inputFile
     
    # DECIMAL RA and DEC VALUES OF HE0435-1223
    ra_qsr = (4.+38./60.+14.9/60./60.)/24*360 
    dec_qsr = -12. - (17./60 +14.4/60./60.)

    npxls = 31
    # INITIALIZE THE CUSTOMIZED FITS FILE MANAGER
    FM = FITSmanager(inputFile)

    # PLOT QSR IMAGE
    FM.plot_image(ra_qsr, dec_qsr, Npx=npxls, out = args.outputFileTag+'_qsr')

    # ESTIMATE THE BACKGROUND AND READ NOISE
    FM.estimate_read_noise()

    # APASS PHOTOMETRY
    APASS_table = ascii.read('../../data/HE0435_LCOGT/APASS_0438_list.csv')
    APASS_rejects = [9, 21, 22] # 9 and 21 seem to be at the edge of the field of view. 22 seems close to the edge, sometimes we don't catch it.
    ZP_mean, ZP_rms, alpha_mean, beta_mean = APASS_zero_points(FM, APASS_table, APASS_rejects, FM.readnoise, display=False)
    print ZP_mean, ZP_rms, alpha_mean, beta_mean

    quadFit(FM, ra_qsr, dec_qsr, ZP_mean, ZP_rms, alpha_mean, beta_mean, npxls)




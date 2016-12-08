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
matplotlib.rcParams['figure.max_open_warning'] = 100

if __name__ == "__main__":

    parser=argparse.ArgumentParser(description='quadPhot routine to calculate crowded field quadruply lensed quasar photometry')
    parser.add_argument("-d",    "--dataDir",        default = '/halo_nobackup/lenssim/romerowo/lcogt_data/Dec_2015_100_hours',              help="directory with FITS file images",type=str)
    parser.add_argument("-i",    "--imageFile",      default = 'coj1m003-kb71-20151214-0083-e91.fits',                                       help="directory with FITS file images",type=str)
    parser.add_argument("-fi",   "--fieldInputs",    default = '/halo_nobackup/lenssim/romerowo/lcolens/python/arw/LCOGT_HE0435_inputs.txt', help="file with field inputs",type=str)
    parser.add_argument("-nmax", "--nmax",           default = None,  help="Maximum n mode for shapelet fitting. Defined in filedInputs. These entries override that for testing and debugging.",type=int)
    parser.add_argument("-mmax", "--mmax",           default = None,  help="Maximum m mode for shapelet fitting. Defined in filedInputs. These entries override that for testing and debugging.",type=int)
    parser.add_argument("-pc",   "--photCal",        default = None,  help="Photometric calibration value can be passed",type=float)
    parser.add_argument("-pcu",  "--photCal_unc",    default = None,  help="Photometric calibration uncertainty",type=float)
    parser.add_argument("-o",    "--outputFileTag",  default = 'out', help="Label for output files. Includes target directory.",type=str)
    parser.add_argument("-p",    "--plots",          default = 3,     help="level of plot production",type=int)
    parser.add_argument("-L",    "--analysis_level", default = 0,     help="analysis levels: PSF=0, Relative Photometry = 1, Crowded Field = 2",type=int)
    parser.add_argument("-pf",   "--PSF_file",       default = None,  help="npz file with PSF star fit results",type=str)
    parser.add_argument("-pcsc", "--PSF_chi_sq_cut", default = None,  help="Maximum allowable PSF chi sq. Defined in filedInputs. These entries override that for testing and debugging.",type=float)
    parser.add_argument("-pmcc", "--PSF_max_chi_cut", default = None,  help="Maximum allowable PSF max(abs(chi)). Defined in filedInputs. These entries override that for testing and debugging.",type=float)


    #KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
    for i, arg in enumerate(sys.argv):
      if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

    args=parser.parse_args()
    imageFile = args.dataDir+'/'+args.imageFile
    print '\n############################################################'
    print 'quadPhot imageFile ', imageFile.split('/')[-1]
    print 'quadPhot fieldInput', args.fieldInputs.split('/')[-1]
    print 'quadPhot outputFile', args.outputFileTag
    print '############################################################\n'
    # print args
    # print 
    print '\nInputs Parameters and Values'
    for arg in vars(args):
        print '\t', arg.ljust(14), '= ', getattr(args, arg)
    field_inputs = eval(open(args.fieldInputs).read())

    npxls = field_inputs['Num_Pixels']

    # Check nmax, mmax inputs
    nmax = field_inputs['nmax']
    mmax = field_inputs['mmax']
    PSF_chi_sq_cut  = field_inputs['PSF_chi_sq_cut']
    PSF_max_chi_cut = field_inputs['PSF_max_chi_cut']
    if(args.nmax!=None and args.mmax!=None):
        nmax = args.nmax
        mmax = args.mmax
    if(args.PSF_chi_sq_cut != None and args.PSF_max_chi_cut != None):
        PSF_chi_sq_cut  = args.PSF_chi_sq_cut
        PSF_max_chi_cut = args.PSF_max_chi_cut


    ''' PSF RUNS '''
    if(args.analysis_level==0):
        # INITIALIZE THE CUSTOMIZED FITS FILE MANAGER
        FM = FITSmanager(imageFile)
        # ESTIMATE THE BACKGROUND AND READ NOISE
        FM.estimate_read_noise(display=args.plots, out = args.outputFileTag+'_readnoise')
        # ESTIMATE THE PSF
        estimate_PSF(FM, field_inputs['PSF_Star_File'], FM.readnoise, nmax, mmax, field_inputs['Num_Pixels'],  display=args.plots, out = args.outputFileTag)
        print 'PSF Estimation Complete'

    ''' Photometric Fitting '''
    if(args.analysis_level==1):
        # GET THE FITS FILE NAME
        f = np.load(args.PSF_file)
        FM = FITSmanager(str(f['imageFile']))
        readnoise = float(f['readnoise'])
        f.close()
        # ESTIMATE THE MEDAIN PSF PARAMETERS
        beta, beta_unc, PSF_parms, PSF_parms_unc, num_good_PSF_fit = Get_Median_PSF_Parameters(args.PSF_file, PSF_chi_sq_cut, PSF_max_chi_cut, nmax, mmax )

        # CHECK S_CCD NORMALIZATION OF PSF PARAMETERS. 
        nm_list = np.array(psf.get_nm(nmax, mmax))
        fn_indices = fn0_indices = np.where(nm_list[:,1]==0)[0]
        #print 'PSF_parms', len(PSF_parms)
        #print 'nm_list', len(nm_list)
        S_CCD = np.sum(PSF_parms[fn0_indices])*(2*np.sqrt(np.pi)*beta)
        print 'normed_parameter S_CCD', S_CCD

        # ESTIMATE THE BRIGHTNESS OF STARS IN THE FIELD
        star_table = ascii.read(field_inputs['Phot_Star_File'])
        #print star_table.keys()
        star_ra  = np.array(star_table['ALPHA_J2000_1'])
        star_dec = np.array(star_table['DELTA_J2000_1'])
	# USE THIS ONE TO CHECK WITH APASS STARS
        #star_ra  = np.array(star_table['radeg'])
        #star_dec = np.array(star_table['decdeg'])
        #star_ra  = np.array(star_table['ALPHA_J2000_2'])
        #star_dec = np.array(star_table['DELTA_J2000_2'])
        #star_ra, star_dec = readStarList2(field_inputs['Phot_Star_File'])
        #print star_ra
        #print star_dec
        fit_success, chi_sq_list, max_chi_list, x0_list, y0_list, bkg_list, S_CCD_list, S_CCD_unc_list = starFit(FM, star_ra, star_dec, beta, PSF_parms, readnoise, nmax, mmax, N_px=npxls, display = args.plots, outputFileTag=args.outputFileTag)
        npz_out = args.outputFileTag + '_Phot.npz'
        np.savez( npz_out,
                starFileName    = field_inputs['Phot_Star_File'],
                imageFile       = FM.fits_file_name,
                outFileTag      = args.outputFileTag,
                beta            = beta, 
                beta_unc        = beta_unc, 
                PSF_parms       = PSF_parms,
                PSF_parms_unc   = PSF_parms_unc, 
                num_good_PSF_fit=num_good_PSF_fit,
                mjd_obs         = float(FM.hdulist[0].header['MJD-OBS']),
                readnoise       = readnoise,
                nmax            = nmax,
                mmax            = mmax,
                filter          = FM.hdulist[0].header['FILTER'],
                S_CCD           = S_CCD_list,
                sig_S_CCD       = S_CCD_unc_list,
                chi_sq          = chi_sq_list, 
                max_chi         = max_chi_list,
                x0              = x0_list,
                y0              = y0_list,
                bkg             = bkg_list,
                fit_success     = fit_success
        )
        print 'Stellar Photometry Complete'

    ''' Lensed Quasar Fitting '''
    if(args.analysis_level==2):
        # GET THE FITS FILE NAME
        f = np.load(args.PSF_file)
        FM = FITSmanager(str(f['imageFile']))
        readnoise = float(f['readnoise'])
        f.close()
        # ESTIMATE THE MEDAIN PSF PARAMETERS
        beta, beta_unc, PSF_parms, PSF_parms_unc, num_good_PSF_fit = Get_Median_PSF_Parameters(args.PSF_file, PSF_chi_sq_cut, PSF_max_chi_cut, nmax, mmax )
        # FIT THE QUAD IMAGE
        # Check that Phot_cal is true (then use expected flux of the galaxy, otherwise, fit for the flux of the lensing galaxy

        # CHECK S_CCD NORMALIZATION OF PSF PARAMETERS. 
        nm_list = np.array(psf.get_nm(nmax, mmax))
        fn_indices = fn0_indices = np.where(nm_list[:,1]==0)[0]
        #print 'PSF_parms', len(PSF_parms)
	    #print 'nm_list', len(nm_list)
        S_CCD = np.sum(PSF_parms[fn0_indices])*(2*np.sqrt(np.pi)*beta)
        print 'normed_parameter S_CCD', S_CCD
        ra_field, dec_field = convertRaAndDec(field_inputs['Field_RA'], field_inputs['Field_Dec'])
        print 'Field RA, DEC', ra_field, dec_field
        popt_wg, pcov_wg, chisq_wg, max_chi_wg = quadFit(FM, float(ra_field), float(dec_field), 
                                                         np.array(field_inputs['Point_Source_RAs']), np.array(field_inputs['Point_Source_Decs']), 
                                                         float(field_inputs['Lens_Gal_RA']), float(field_inputs['Lens_Gal_Dec']), 
                                                         readnoise, beta, PSF_parms, nmax, mmax, 
                                                         npxls, galFit=True, display=args.plots,  
                                                         outputFileTag=args.outputFileTag, emcee_level = 0)
    #########################################################################################

    if(args.analysis_level>2):

        print field_inputs
        # DECIMAL RA and DEC VALUES OF HE0435-1223
        ra_field, dec_field = convertRaAndDec(field_inputs['Field_RA'], field_inputs['Field_Dec'])
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

        popt_ng, pcov_ng, chisq_ng, max_chi_ng = quadFit(FM, ra_qsr, dec_qsr, ra_images, dec_images, ra_lensgal, dec_lensgal, beta_mean, shapelet_coeffs, npxls, galFit=False, display=args.plots,  outputFileTag=args.outputFileTag, emcee_level = args.emcee_level)


        popt_wg, pcov_wg, chisq_wg, max_chi_wg = quadFit(FM, ra_qsr, dec_qsr, ra_images, dec_images, ra_lensgal, dec_lensgal, beta_mean, shapelet_coeffs, npxls, galFit=True, display=args.plots,  outputFileTag=args.outputFileTag, emcee_level = args.emcee_level)

        #star_ra, star_dec = readStarList(args.star_list)
        #star_index_list, star_chi_sq, star_max_chi, star_S_CCD, star_S_CCD_unc = starFit(FM, star_ra, star_dec, beta_mean, shapelet_coeffs, N_px=npxls, display = args.plots, outputFileTag=args.outputFileTag)

        #exit()

        filter   = FM.hdulist[0].header['FILTER']
        mjd_obs  = float(FM.hdulist[0].header['MJD-OBS'])
        npz_out  = args.outputFileTag + '_results.npz'

        np.savez(npz_out,
	    imageFile       = imageFile,
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
        qsr_wg_chisq    = chisq_wg, 
        qsr_wg_max_chi  = max_chi_wg,
        qsr_ng_parms    = popt_ng, 
        qsr_ng_covar    = pcov_ng,
        qsr_ng_chisq    = chisq_ng, 
        qsr_ng_max_chi  = max_chi_ng
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



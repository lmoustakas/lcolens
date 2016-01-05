import matplotlib
matplotlib.use('Agg')
import numpy as np
from pylab import *
from astropy.io import ascii
import glob
import os
import AnalysisLibrary as AL
from time import clock

rcParams['font.size']=24
rcParams['legend.fontsize']=18
rcParams['figure.facecolor']='white'
#fnames = glob.glob('/nisushome/romerowo/lcolens_20150605/python/arw/npzfiles/image_*_results.npz') # no priors of seeing imposed
#fnames = glob.glob('/disk4/romerowo/lcolens_outputs/20150916/npzfiles/image_*_results.npz') # priors of seeing imposed
#fnames = glob.glob('/disk4/romerowo/lcolens_outputs/20150917/npzfiles/image_*_results.npz') # added pixelization and robust ZP (botched)
#fnames = glob.glob('/disk4/romerowo/lcolens_outputs/20150918/npzfiles/image_*_results.npz') # robust ZP, no subsample-integrated pixelization
#fnames = glob.glob('/disk4/romerowo/lcolens_outputs/20151013/npzfiles/image_*_results.npz') # robust ZP, no subsample-integrated pixelization
fnames = glob.glob('/disk4/romerowo/lcolens_outputs/20151020/npzfiles/image_*_results.npz') # robust ZP, added subsample-integrated pixelization
#fnames = glob.glob('/disk4/romerowo/lcolens_outputs/20151028/npzfiles/image_*_results.npz') # robust ZP, fixed bug and removed subsample-integrated pixelization

#fnames = glob.glob('/nisushome/romerowo/lcolens_20150605/python/arw/npzfiles/image_21*_results.npz')


fnames = sorted(fnames)
for f in fnames:
	print f

def read_emcee_chains(samples):
	nwalkers = 100
	niterations=10000
	def condition_chain(k):
		return np.ravel(samples[:,k].reshape(nwalkers, niterations)[:,niterations/2:])
		#return np.ravel(samples[:,k].reshape(nwalkers, niterations)[:,niterations/2:])

	x0Ch    = condition_chain(0)
	y0Ch    = condition_chain(1) 
	amp1Ch  = condition_chain(2)
	amp2Ch  = condition_chain(3)
	amp3Ch  = condition_chain(4)
	amp4Ch  = condition_chain(5)
	alphaCh = condition_chain(6)
	betaCh  = condition_chain(7)
	nbkgCh  = condition_chain(8)
	return x0Ch, y0Ch, amp1Ch, amp2Ch, amp3Ch, amp4Ch, alphaCh, betaCh, nbkgCh


mjd_obs=[]
m1 = []
me1 = []
m2 = []
me2 = []
m3 = []
me3 = []
m4 = []
me4 = []
maxChi = []
chiSq = []
readnoise=[]
#readnoise_error=[]
site=[]
ZP_mean=[]
ZP_rms=[]
ZP_wrms=[]
ZP_std=[]
filter=[]
ZP_flx = []
alpha=[]
beta=[]
fwhmEmcee=[]
fwhmEmceeUpperError=[]
fwhmEmceeLowerError=[]
amp1Emcee=[]
amp1EmceeUpperError=[]
amp1EmceeLowerError=[]
amp2Emcee=[]
amp2EmceeUpperError=[]
amp2EmceeLowerError=[]
amp3Emcee=[]
amp3EmceeUpperError=[]
amp3EmceeLowerError=[]
amp4Emcee=[]
amp4EmceeUpperError=[]
amp4EmceeLowerError=[]
nbkgEmcee=[]
nbkgEmceeUpperError=[]
nbkgEmceeLowerError=[]
LC1Emcee=[]
LC1EmceeUpperError=[]
LC1EmceeLowerError=[]
LC2Emcee=[]
LC2EmceeUpperError=[]
LC2EmceeLowerError=[]
LC3Emcee=[]
LC3EmceeUpperError=[]
LC3EmceeLowerError=[]
LC4Emcee=[]
LC4EmceeUpperError=[]
LC4EmceeLowerError=[]
airmass = []
seeing_fwhm = []
input_files = []
det_cov = []
pxscl = []
fnm_list=[]
outFileTags=[]
theta_list = []


#rand_vals = arange(np.random.randint(0,len(fnames))-5,5)
#rand_vals = np.random.randint(0,len(fnames),5)
#rand_vals = sort(rand_vals)
#for k in [30, 80,90,100, 150, 210]:
#for k in rand_vals:
#for k in range(0,2):
for k in range(0,len(fnames)):
    start_time = clock()
    print '\n\n############################################\n\n'
    print 'k=',k, fnames[k]
    if(os.path.exists(fnames[k].replace('image', 'bad_image'))): continue # bad covariance reconstruction

    results = np.load(fnames[k])
    #print '\t',str(results['inputFile'])
    #if('/coj' not in str(results['inputFile'])): continue

    mjd_obs.append(results['mjd_obs'])
    m1.append(results['m1'])
    me1.append(results['me1'])
    m2.append(results['m2'])
    me2.append(results['me2'])
    m3.append(results['m3'])
    me3.append(results['me3'])
    m4.append(results['m4'])
    me4.append(results['me4'])
    maxChi.append(results['maxChi'])
    chiSq.append(results['chiSq'])
    readnoise.append(results['readnoise'])
    #readnoise_error.append(results['readnoise_error'])
    site.append(str(results['inputFile']).split('/')[-1][0:3])
    ZP_mean.append(results['ZP_mean'])
    ZP_rms.append(results['ZP_rms'])
    ZP_wrms.append(results['ZP_wrms'])
    #ZP_std.append(results['ZP_std'])
    filter.append(results['filter'])
    alpha.append(results['APASS_alpha'])
    beta.append(results['APASS_beta'])
    outFileTags.append(results['outFileTag'])
    #if results['chiSq'] > 2. or results['maxChi'] > 5.:
    if(k%20==0): print 'outFnm\t\tmjd_obs\t\tchiSq\tmaxChi\tinputFile\t\t\t\tm1\tZP_mean\tfilter'
    print '%s\t%1.3f\t%1.2f\t%1.2f\t%s\t%1.2f\t%1.2f\t%s'%(results['outFileTag'], results['mjd_obs'], results['chiSq'], results['maxChi'], str(results['inputFile']).split('/')[-1], results['m1'], results['ZP_mean'], results['filter'])
    input_files.append(str(results['inputFile']).split('/')[-1])
    if k==0:
		print '\t',results.files
    fnm = str(results['inputFile'])
    fnm_list.append(fnm)
    #fnm = '/nisushome'+fnm
    print '\t',fnm
    print '\tloading chains: %1.2f s'%(clock()-start_time)
    chain_fnm = fnames[k].replace('results', 'chains')
    fchain = np.load(chain_fnm)
    print '\testimating means and uncertainties chains: %1.2f s'%(clock()-start_time)
    samples = fchain['arr_0']
    x0Ch, y0Ch, amp1Ch, amp2Ch, amp3Ch, amp4Ch, alphaCh, betaCh, nbkgCh = read_emcee_chains(samples)
    fwhmCh = alphaCh*2.*np.sqrt(2.**(1/betaCh)-1.)
    parms = lambda v: (v[1], v[2]-v[1], v[1]-v[0]), (np.percentile(fwhmCh, [16, 50, 84],axis=0))
    fwhmEmcee.append(parms[1][1])
    fwhmEmceeUpperError.append(parms[1][1]-parms[1][0])
    fwhmEmceeLowerError.append(parms[1][2]-parms[1][1])

    parms = lambda v: (v[1], v[2]-v[1], v[1]-v[0]), (np.percentile(amp1Ch, [16, 50, 84],axis=0))
    amp1Emcee.append(parms[1][1])
    amp1EmceeUpperError.append(parms[1][1]-parms[1][0])
    amp1EmceeLowerError.append(parms[1][2]-parms[1][1])
    parms = lambda v: (v[1], v[2]-v[1], v[1]-v[0]), (np.percentile(amp2Ch, [16, 50, 84],axis=0))
    amp2Emcee.append(parms[1][1])
    amp2EmceeUpperError.append(parms[1][1]-parms[1][0])
    amp2EmceeLowerError.append(parms[1][2]-parms[1][1])
    parms = lambda v: (v[1], v[2]-v[1], v[1]-v[0]), (np.percentile(amp3Ch, [16, 50, 84],axis=0))
    amp3Emcee.append(parms[1][1])
    amp3EmceeUpperError.append(parms[1][1]-parms[1][0])
    amp3EmceeLowerError.append(parms[1][2]-parms[1][1])
    parms = lambda v: (v[1], v[2]-v[1], v[1]-v[0]), (np.percentile(amp4Ch, [16, 50, 84],axis=0))
    amp4Emcee.append(parms[1][1])
    amp4EmceeUpperError.append(parms[1][1]-parms[1][0])
    amp4EmceeLowerError.append(parms[1][2]-parms[1][1])
    parms = lambda v: (v[1], v[2]-v[1], v[1]-v[0]), (np.percentile(nbkgCh, [16, 50, 84],axis=0))
    nbkgEmcee.append(parms[1][1])
    nbkgEmceeUpperError.append(parms[1][1]-parms[1][0])
    nbkgEmceeLowerError.append(parms[1][2]-parms[1][1])


    dirc = '/nisushome/data2/romerowo/lcogt_data/he045-1223_wcs_corrected/'
    #print '%s'%(fnm)
    #exit()
    print '\topening FITS file: %1.2f s'%(clock()-start_time)
    FM = AL.FITSmanager(fnm)
    airmass.append(FM.hdulist[0].header['AIRMASS'])


    light_distrib1=[]
    light_distrib2=[]
    light_distrib3=[]
    light_distrib4=[]
    xg,yg = np.mgrid[:31, :31]
        
    print 'random resample of light estimated: %1.2f s'%(clock()-start_time)
    for j in np.random.randint(0,len(x0Ch),1000):
		#img = quad_image_model(FM, x0Ch[j], y0Ch[j], amp1Ch[j], amp2Ch[j], amp3Ch[j], amp4Ch[j], alphaCh[j], betaCh[j], nbkgCh[j], N_pix=31, flip=False)
		theta1 = [amp1Ch[j], alphaCh[j], betaCh[j], 15, 15, nbkgCh[j]]
		theta2 = [amp2Ch[j], alphaCh[j], betaCh[j], 15, 15, nbkgCh[j]]
		theta3 = [amp3Ch[j], alphaCh[j], betaCh[j], 15, 15, nbkgCh[j]]
		theta4 = [amp4Ch[j], alphaCh[j], betaCh[j], 15, 15, nbkgCh[j]]
		img1 = AL.twoD_Moffat((xg, yg), *theta1).reshape(31,31)
		img2 = AL.twoD_Moffat((xg, yg), *theta2).reshape(31,31)
		img3 = AL.twoD_Moffat((xg, yg), *theta3).reshape(31,31)
		img4 = AL.twoD_Moffat((xg, yg), *theta4).reshape(31,31)
		intg1 = np.sum(img1 - nbkgCh[j])
		intg2 = np.sum(img2 - nbkgCh[j])
		intg3 = np.sum(img3 - nbkgCh[j])
		intg4 = np.sum(img4 - nbkgCh[j])
		light_distrib1.append(intg1)
		light_distrib2.append(intg2)
		light_distrib3.append(intg3)
		light_distrib4.append(intg4)
		#print j, len(x0Ch), '%1.2e\t%1.2e\t%1.2e\t%1.2e'%(intg1,intg2,intg3,intg4)
    '''
    subplot(221)
    hist(light_distrib1)
    subplot(222)
    hist(light_distrib2)
    subplot(223)
    hist(light_distrib3)
    subplot(224)
    hist(light_distrib4)
    '''
    print 'estimating magnitudes: %1.2f s'%(clock()-start_time)
    parms = lambda v: (v[1], v[2]-v[1], v[1]-v[0]), (np.percentile(light_distrib1, [16, 50, 84],axis=0))
    #print parms
    LC1Emcee.append(parms[1][1])
    LC1EmceeUpperError.append(parms[1][1]-parms[1][0])
    LC1EmceeLowerError.append(parms[1][2]-parms[1][1])
    parms = lambda v: (v[1], v[2]-v[1], v[1]-v[0]), (np.percentile(light_distrib2, [16, 50, 84],axis=0))
    #print parms
    LC2Emcee.append(parms[1][1])
    LC2EmceeUpperError.append(parms[1][1]-parms[1][0])
    LC2EmceeLowerError.append(parms[1][2]-parms[1][1])
    parms = lambda v: (v[1], v[2]-v[1], v[1]-v[0]), (np.percentile(light_distrib3, [16, 50, 84],axis=0))
    #print parms
    LC3Emcee.append(parms[1][1])
    LC3EmceeUpperError.append(parms[1][1]-parms[1][0])
    LC3EmceeLowerError.append(parms[1][2]-parms[1][1])
    parms = lambda v: (v[1], v[2]-v[1], v[1]-v[0]), (np.percentile(light_distrib4, [16, 50, 84],axis=0))
    #print parms
    LC4Emcee.append(parms[1][1])
    LC4EmceeUpperError.append(parms[1][1]-parms[1][0])
    LC4EmceeLowerError.append(parms[1][2]-parms[1][1])
    #show()

    print '\tLCEmcee', LC1Emcee[-1], LC2Emcee[-1], LC3Emcee[-1], LC4Emcee[-1]
    # DECIMAL RA and DEC VALUES OF HE0435-1223
    ra_qsr = (4.+38./60.+14.9/60./60.)/24*360 
    dec_qsr = -12. - (17./60 +14.4/60./60.)
    #FM = AL.FITSmanager(fnm) 

    ZP_flux       = np.array(AL.magnitude2flux(np.array(ZP_mean)))
    #ZP_flux_std   = np.array(AL.magErr2fluxErr(np.array(ZP_mean), np.array(ZP_std)))
    ZP_flux_wrms  = np.array(AL.magErr2fluxErr(np.array(ZP_mean), np.array(ZP_wrms)))
    flux1 = np.array(ZP_flux*LC1Emcee)
    flux1UpperError=np.array(flux1*np.sqrt(np.array(LC1EmceeUpperError)**2/np.array(LC1Emcee)**2 + ZP_flux_wrms**2/ZP_flux**2))
    flux1LowerError=np.array(flux1*np.sqrt(np.array(LC1EmceeLowerError)**2/np.array(LC1Emcee)**2 + ZP_flux_wrms**2/ZP_flux**2))
    flux2 = np.array(ZP_flux*LC2Emcee)
    flux2UpperError=np.array(flux2*np.sqrt(np.array(LC2EmceeUpperError)**2/np.array(LC2Emcee)**2 + ZP_flux_wrms**2/ZP_flux**2))
    flux2LowerError=np.array(flux2*np.sqrt(np.array(LC2EmceeLowerError)**2/np.array(LC2Emcee)**2 + ZP_flux_wrms**2/ZP_flux**2))
    flux3 = np.array(ZP_flux*LC3Emcee)
    flux3UpperError=np.array(flux3*np.sqrt(np.array(LC3EmceeUpperError)**2/np.array(LC3Emcee)**2 + ZP_flux_wrms**2/ZP_flux**2))
    flux3LowerError=np.array(flux3*np.sqrt(np.array(LC3EmceeLowerError)**2/np.array(LC3Emcee)**2 + ZP_flux_wrms**2/ZP_flux**2))
    flux4 = np.array(ZP_flux*LC4Emcee)
    flux4UpperError=np.array(flux4*np.sqrt(np.array(LC4EmceeUpperError)**2/np.array(LC4Emcee)**2 + ZP_flux_wrms**2/ZP_flux**2))
    flux4LowerError=np.array(flux4*np.sqrt(np.array(LC4EmceeLowerError)**2/np.array(LC4Emcee)**2 + ZP_flux_wrms**2/ZP_flux**2))

    mag1 = AL.flux2magnitude(flux1)
    me1Upper = AL.fluxErr2magErr(flux1, flux1UpperError)
    me1Lower = AL.fluxErr2magErr(flux1, flux1LowerError)
    mag2 = AL.flux2magnitude(flux2)
    me2Upper = AL.fluxErr2magErr(flux2, flux2UpperError)
    me2Lower = AL.fluxErr2magErr(flux2, flux2LowerError)
    mag3 = AL.flux2magnitude(flux3)
    me3Upper = AL.fluxErr2magErr(flux3, flux3UpperError)
    me3Lower = AL.fluxErr2magErr(flux3, flux3LowerError)
    mag4 = AL.flux2magnitude(flux4)
    me4Upper = AL.fluxErr2magErr(flux4, flux4UpperError)
    me4Lower = AL.fluxErr2magErr(flux4, flux4LowerError)

    mag1_array = AL.flux2magnitude(ZP_flux[-1]*np.array(light_distrib1))
    mag2_array = AL.flux2magnitude(ZP_flux[-1]*np.array(light_distrib2))
    mag3_array = AL.flux2magnitude(ZP_flux[-1]*np.array(light_distrib3))
    mag4_array = AL.flux2magnitude(ZP_flux[-1]*np.array(light_distrib4))
    # estimate covariance matrix
    mu1 = np.mean(mag1_array)
    mu2 = np.mean(mag2_array)
    mu3 = np.mean(mag3_array)
    mu4 = np.mean(mag4_array)
    print 'estimating covariance matrix: %1.2f s'%(clock()-start_time)
    print ZP_flux[-1]
    for ii in range(0,len(light_distrib2)):
        if light_distrib2[ii]<=0.:
            print 'negative value!', ii, light_distrib2[ii]
    print np.sum(np.isnan(ZP_flux[-1]*np.array(light_distrib1))), len(ZP_flux[-1]*np.array(light_distrib1))
    print np.sum(np.isnan(ZP_flux[-1]*np.array(light_distrib2))), len(ZP_flux[-1]*np.array(light_distrib2))
    print np.sum(np.isnan(ZP_flux[-1]*np.array(light_distrib3))), len(ZP_flux[-1]*np.array(light_distrib3))
    print np.sum(np.isnan(ZP_flux[-1]*np.array(light_distrib4))), len(ZP_flux[-1]*np.array(light_distrib4))
    print ''
    print np.sum(np.isnan(mag1_array)), len(mag1_array)
    print np.sum(np.isnan(mag2_array)), len(mag2_array)
    print np.sum(np.isnan(mag3_array)), len(mag3_array)
    print np.sum(np.isnan(mag4_array)), len(mag4_array)
    cov     = np.cov([mag1_array, mag2_array, mag3_array, mag4_array])
    cov_det = np.linalg.det(cov)
    print '\t','***** COV DET **** %1.2e'%cov_det
    icov = np.infty
    icov_det = np.infty

    #cov_mat_det = np.linalg.det(cov_mat)
    print cov
    print cov_det
    try: 
        	icov = np.linalg.inv(cov)
        	icov_det = np.linalg.det(icov)
    except:
		print '\t !!!!!!! Could not invert covariance matrix!'
		exit()

    print '\t','***** ICOV DET ****%1.2e'%icov_det
    det_cov.append(cov_det)
    print '\tmag1_array_stats:',np.mean(mag1_array), np.sqrt(np.mean((mag1_array - mu1)**2)), np.sqrt(cov_det)**(1./4.)
    print '\tmag2_array_stats:',np.mean(mag2_array), np.sqrt(np.mean((mag2_array - mu1)**2)), np.sqrt(cov_det)**(1./4.)
    print '\tmag3_array_stats:',np.mean(mag3_array), np.sqrt(np.mean((mag3_array - mu1)**2)), np.sqrt(cov_det)**(1./4.)
    print '\tmag4_array_stats:',np.mean(mag4_array), np.sqrt(np.mean((mag4_array - mu1)**2)), np.sqrt(cov_det)**(1./4.)
    #exit()

    mjd0 = np.floor(np.min(np.array(mjd_obs)))
    seeing_fwhm.append(np.mean(alphaCh)*2.*np.sqrt(2.**(1/np.mean(betaCh))-1.))
    pxscl.append(float(FM.hdulist[0].header['PIXSCALE']))

    if(icov_det==np.infty):
		continue
    print '\tZP_flux, ZP_mean',AL.magnitude2flux(results['ZP_mean']),results['ZP_mean']
    print '\tx0,y0',np.mean(x0Ch), np.mean(y0Ch)
    ZP_flx.append(AL.magnitude2flux(results['ZP_mean']))
    theta = [np.mean(x0Ch),np.mean(y0Ch),np.mean(amp1Ch),np.mean(amp2Ch),np.mean(amp3Ch),np.mean(amp4Ch), np.mean(alphaCh), np.mean(betaCh), np.mean(nbkgCh)]
    theta_list.append(theta)
    print 'plotting figures: %1.2f s'%(clock()-start_time)

    FM.hdulist.close()
    results.close()
    del FM
    print 'loop iteration finished: %1.2f s'%(clock()-start_time)

mjd_obs   = array(mjd_obs)
mag1 = np.array(mag1)
mag2 = np.array(mag2)
mag3 = np.array(mag3)
mag4 = np.array(mag4)
mag1_err = (me1Upper+me1Lower)/2.
mag2_err = (me2Upper+me2Lower)/2.
mag3_err = (me3Upper+me3Lower)/2.
mag4_err = (me4Upper+me4Lower)/2.

###################################
# SAVE FULL RESOLUTION FIGURES ####
###################################

seeing_arcsec = np.array(seeing_fwhm)*np.array(pxscl)
cov_uncertainty = np.power(np.sqrt(det_cov),1./4.)
np.savez('emceeResults.npz', mjd_obs=mjd_obs, airmass=airmass, seeing_arcsec=seeing_arcsec, cov_uncertainty=cov_uncertainty, mag1=mag1, mag2=mag2, mag3=mag3, mag4=mag4, me1Lower=me1Lower, me2Lower=me2Lower, me3Lower=me3Lower, me4Lower=me4Lower, me1Upper=me1Upper, me2Upper=me2Upper, me3Upper=me3Upper, me4Upper=me4Upper, input_files=input_files)
npzfile = np.load('emceeResults.npz')
print npzfile.files

figure(1)
ax = subplot(111)
for nn in range(0,len(mjd_obs)):
	if 'lsc' in input_files[nn]:
		plt.plot([mjd_obs[nn] - mjd0], [airmass[nn]], 'bs')
	if 'coj' in input_files[nn]:
		plt.plot([mjd_obs[nn] - mjd0], [airmass[nn]], 'go')
	if 'cpt' in input_files[nn]:
		plt.plot([mjd_obs[nn] - mjd0], [airmass[nn]], 'r^')
#if 'lsc' in input_files[k]:
#	plt.plot([mjd_obs[k] - mjd0], [airmass[k]], 's', ms=15)
#if 'coj' in input_files[k]:
#	plt.plot([mjd_obs[k] - mjd0], [airmass[k]], 'go', ms=15)
#if 'cpt' in input_files[k]:
#	plt.plot([mjd_obs[k] - mjd0], [airmass[k]], 'r^', ms=15)
plt.plot([-1.], [-1.], 'bs', ms=8, label='Chile')
plt.plot([-1.], [-1.], 'go', ms=8, label='Australia')
plt.plot([-1.], [-1.], 'r^', ms=8, label='S. Africa')
plt.legend(loc=(0.17,0.), numpoints=1, fontsize=24, frameon=False, borderaxespad=-0.05, handletextpad=-0.5, columnspacing=-0.5, labelspacing=0., borderpad=-0.05)
#plot([mjd_obs[k] - mjd0, mjd_obs[k] - mjd0], [2.3, 1.0], 'k--')
plt.xlim(0.,2.3)
plt.ylim(2.3,1.0)
y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.yaxis.set_major_formatter(y_formatter)
ylabel('Airmass')
plt.savefig('airmass.pdf')


figure(2)
for nn in range(0,len(mjd_obs)):
	if 'lsc' in input_files[nn]:
		plt.plot([mjd_obs[nn] - mjd0], [seeing_fwhm[nn]*pxscl[nn]], 'bs')
	if 'cpt' in input_files[nn]:
		plt.plot([mjd_obs[nn] - mjd0], [seeing_fwhm[nn]*pxscl[nn]], 'r^')
	if 'coj' in input_files[nn]:
		plt.plot([mjd_obs[nn] - mjd0], [seeing_fwhm[nn]*pxscl[nn]], 'go')
#if 'lsc' in input_files[k]:
#	plt.plot([mjd_obs[k] - mjd0], [seeing_fwhm[k]*pxscl[k]], 'bs', ms=15)
#if 'cpt' in input_files[k]:
#	plt.plot([mjd_obs[k] - mjd0], [seeing_fwhm[k]*pxscl[k]], 'r^', ms=15)
#if 'coj' in input_files[k]:
#	plt.plot([mjd_obs[k] - mjd0], [seeing_fwhm[k]*pxscl[k]], 'go', ms=15)
# plt.plot(mjd_obs - mjd0, seeing_fwhm, 'ro')
#plot([mjd_obs[k] - mjd0,mjd_obs[k] - mjd0], [0.,2.75], 'k--')
plt.xlim(0.,2.3)
plt.ylim(0.,2.75)
y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.yaxis.set_major_formatter(y_formatter)
ylabel('Seeing, arcseconds')
plt.savefig('seeing.pdf')

figure(3, figsize=(8,10))
ax = subplot(111)
ax.set_yscale('log')
for nn in range(0,len(mjd_obs)):
	if 'lsc' in input_files[nn]:
		plt.plot([mjd_obs[nn] - mjd0], [np.power(np.sqrt(det_cov[nn]),1./4.)], 'bs')
	if 'cpt' in input_files[nn]:
		plt.plot([mjd_obs[nn] - mjd0], [np.power(np.sqrt(det_cov[nn]),1./4.)], 'r^')
	if 'coj' in input_files[nn]:
		plt.plot([mjd_obs[nn] - mjd0], [np.power(np.sqrt(det_cov[nn]),1./4.)], 'go')
#if 'lsc' in input_files[k]:
#	plt.plot([mjd_obs[k] - mjd0], [np.power(np.sqrt(det_cov[k]),1./4.)], 'bs', ms=15)
#if 'cpt' in input_files[k]:
#	plt.plot([mjd_obs[k] - mjd0], [np.power(np.sqrt(det_cov[k]),1./4.)], 'r^', ms=15)
#if 'coj' in input_files[k]:
#	plt.plot([mjd_obs[k] - mjd0], [np.power(np.sqrt(det_cov[k]),1./4.)], 'go', ms=15)
#plt.plot(mjd_obs - mjd0, np.log10(np.array(det_cov)), 'ro')
#plot([mjd_obs[k] - mjd0,mjd_obs[k] - mjd0], [0.01,0.5], 'k--')
plt.xlim(0.,2.3)
plt.ylim(0.01,0.5)
#y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
#ax.yaxis.set_major_formatter(y_formatter)
ylabel('$(\sqrt{|\Sigma|})^{1/4}$, mag')
xlabel('Days Since MJD %1.0f'%mjd0)
yticks([0.003, 0.01, 0.03, 0.10, 0.3])
title('Magnitude Uncertainty')
grid(True)
subplots_adjust(left=0.2, bottom=0.1)
plt.savefig('covariance.pdf')

figure(4, figsize=(8,10))
errorbar(mjd_obs - mjd0, mag1-1.0, yerr=(me1Lower, me1Upper), fmt=',', color=[0.,0.,0.], elinewidth=2)
#errorbar([mjd_obs[k] - mjd0], [mag1[k]-1.0], yerr=([me1Lower[k]], [me1Upper[k]]), fmt=',', color=[0.,0.,0.],elinewidth=10)  #, label='M1 - 1.0'

errorbar(mjd_obs - mjd0, mag3-0.5, yerr=(me3Lower, me3Upper), fmt=',', color=[0.,0.,0.8], elinewidth=2)
#errorbar([mjd_obs[k] - mjd0], [mag3[k]-0.5], yerr=([me3Lower[k]], [me3Upper[k]]), fmt=',', color=[0.,0.,0.8],elinewidth=10) #, label='M2 - 0.5'

errorbar(mjd_obs - mjd0, mag2, yerr=(me2Lower, me2Upper), fmt=',', color=[0.7,0.,0.], elinewidth=2)
#errorbar([mjd_obs[k] - mjd0], [mag2[k]], yerr=([me2Lower[k]], [me2Upper[k]]), fmt=',', color=[0.7,0.,0.],elinewidth=10)     #, label='S1 + 0.0'

errorbar(mjd_obs - mjd0, mag4+0.5, yerr=(me4Lower, me4Upper), fmt=',', color='g', elinewidth=2)
#errorbar([mjd_obs[k] - mjd0], [mag4[k]+0.5], yerr=([me4Lower[k]], [me4Upper[k]]), fmt=',', color='g',elinewidth=10) #, label='S2 + 0.5'

plt.errorbar([-1.], [-1.], yerr=[1.], fmt='k,', elinewidth=10, label='M1 - 1.0')
plt.errorbar([-1.], [-1.], yerr=[1.], fmt='b,',elinewidth=10,  label='M2 - 0.5')
plt.errorbar([-1.], [-1.], yerr=[1.], fmt='r,',elinewidth=10,  label='S1 + 0.0')
plt.errorbar([-1.], [-1.], yerr=[1.], fmt='g,',elinewidth=10,  label='S2 + 0.5')


#plot([mjd_obs[k] - mjd0,mjd_obs[k] - mjd0], [20.3,16.], 'k--')
plt.xlim(0.,2.3)
plt.ylim(20.3,16.)
y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.yaxis.set_major_formatter(y_formatter)
xlabel('Days Since MJD %1.0f'%mjd0)
ylabel('Image Magnitudes')
legend(loc=2, numpoints=1)
suptitle('HE0435-1223 \n LCOGT Observations', fontsize=36)
grid(True)
subplots_adjust(left=0.15, top=0.85)
plt.savefig('light_curves.pdf')

#suptitle(fnm.split('/')[-1])
#exit()


#for k in range(0,10):
for k in range(0,len(fnames)):
    start_time = clock()
    print '\n\n############################################\n\n'
    print 'k=',k, fnames[k]
    if(os.path.exists(fnames[k].replace('image', 'bad_image'))): continue # bad covariance reconstruction

    print '\topening FITS file: %1.2f s'%(clock()-start_time)
    FM = AL.FITSmanager(fnm_list[k])

    print '\tLCEmcee', LC1Emcee[-1], LC2Emcee[-1], LC3Emcee[-1], LC4Emcee[-1]
    # DECIMAL RA and DEC VALUES OF HE0435-1223
    ra_qsr = (4.+38./60.+14.9/60./60.)/24*360 
    dec_qsr = -12. - (17./60 +14.4/60./60.)
    #FM = AL.FITSmanager(fnm) 


    mjd0 = np.floor(np.min(np.array(mjd_obs)))
    seeing_fwhm.append(np.mean(alphaCh)*2.*np.sqrt(2.**(1/np.mean(betaCh))-1.))
    pxscl.append(float(FM.hdulist[0].header['PIXSCALE']))

    figure(figsize=(19.,24))
    ax1=plt.subplot(322)
    ax2=plt.subplot(324)
    ax3=plt.subplot(326)
    print 'plotting figures: %1.2f s'%(clock()-start_time)
    # FROM CASTLES
    ra_images=np.array([0.,-1.476,-2.467,-0.939])
    dec_images=np.array([0.,0.553, -0.603, -1.614])
    ra_images -= np.mean(ra_images)
    dec_images -= np.mean(dec_images)

    #FM.plot_image_movie(ra_qsr, dec_qsr, ax1, ax2, ax3, ZP_flx[k], theta_list[k], Npx=31)
    FM.plot_image_movie(ra_qsr, dec_qsr, ra_images, dec_images, ax1, ax2, ax3, ZP_flx[k], theta_list[k], Npx=31)
    #FM.plot_image_movie(ra_qsr, dec_qsr, ra_images, dec_images, ax1, ax2, ax3, ZP_flx, theta, Npx=31)
    plt.subplots_adjust(left=0.09)

    #print '\tPIXSCALE',pxscl
    ax=subplot(6,2,1)
    for nn in range(0,len(mjd_obs)):
		if 'lsc' in input_files[nn]:
			plt.plot([mjd_obs[nn] - mjd0], [airmass[nn]], 'bs')
		if 'coj' in input_files[nn]:
			plt.plot([mjd_obs[nn] - mjd0], [airmass[nn]], 'go')
		if 'cpt' in input_files[nn]:
			plt.plot([mjd_obs[nn] - mjd0], [airmass[nn]], 'r^')
    if 'lsc' in input_files[k]:
		plt.plot([mjd_obs[k] - mjd0], [airmass[k]], 's', ms=15)
    if 'coj' in input_files[k]:
		plt.plot([mjd_obs[k] - mjd0], [airmass[k]], 'go', ms=15)
    if 'cpt' in input_files[k]:
		plt.plot([mjd_obs[k] - mjd0], [airmass[k]], 'r^', ms=15)
    plt.plot([-1.], [-1.], 'bs', ms=8, label='Chile')
    plt.plot([-1.], [-1.], 'go', ms=8, label='Australia')
    plt.plot([-1.], [-1.], 'r^', ms=8, label='S. Africa')
    plt.legend(loc=(0.17,0.), numpoints=1, fontsize=24, frameon=False, borderaxespad=-0.05, handletextpad=-0.5, columnspacing=-0.5, labelspacing=0., borderpad=-0.05)
    plot([mjd_obs[k] - mjd0, mjd_obs[k] - mjd0], [2.3, 1.0], 'k--')
    plt.xlim(0.,2.3)
    plt.ylim(2.3,1.0)
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
    ylabel('Airmass')

    ax=subplot(6,2,3)
    for nn in range(0,len(mjd_obs)):
		if 'lsc' in input_files[nn]:
			plt.plot([mjd_obs[nn] - mjd0], [seeing_fwhm[nn]*pxscl[nn]], 'bs')
		if 'cpt' in input_files[nn]:
			plt.plot([mjd_obs[nn] - mjd0], [seeing_fwhm[nn]*pxscl[nn]], 'r^')
		if 'coj' in input_files[nn]:
			plt.plot([mjd_obs[nn] - mjd0], [seeing_fwhm[nn]*pxscl[nn]], 'go')
    if 'lsc' in input_files[k]:
		plt.plot([mjd_obs[k] - mjd0], [seeing_fwhm[k]*pxscl[k]], 'bs', ms=15)
    if 'cpt' in input_files[k]:
		plt.plot([mjd_obs[k] - mjd0], [seeing_fwhm[k]*pxscl[k]], 'r^', ms=15)
    if 'coj' in input_files[k]:
		plt.plot([mjd_obs[k] - mjd0], [seeing_fwhm[k]*pxscl[k]], 'go', ms=15)
    # plt.plot(mjd_obs - mjd0, seeing_fwhm, 'ro')
    plot([mjd_obs[k] - mjd0,mjd_obs[k] - mjd0], [0.,2.75], 'k--')
    plt.xlim(0.,2.3)
    plt.ylim(0.,2.75)
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
    ylabel('Seeing, arcseconds')

    ax=subplot(6,2,5)
    ax.set_yscale('log')
    for nn in range(0,len(mjd_obs)):
		if 'lsc' in input_files[nn]:
			plt.plot([mjd_obs[nn] - mjd0], [np.power(np.sqrt(det_cov[nn]),1./4.)], 'bs')
		if 'cpt' in input_files[nn]:
			plt.plot([mjd_obs[nn] - mjd0], [np.power(np.sqrt(det_cov[nn]),1./4.)], 'r^')
		if 'coj' in input_files[nn]:
			plt.plot([mjd_obs[nn] - mjd0], [np.power(np.sqrt(det_cov[nn]),1./4.)], 'go')
    if 'lsc' in input_files[k]:
		plt.plot([mjd_obs[k] - mjd0], [np.power(np.sqrt(det_cov[k]),1./4.)], 'bs', ms=15)
    if 'cpt' in input_files[k]:
		plt.plot([mjd_obs[k] - mjd0], [np.power(np.sqrt(det_cov[k]),1./4.)], 'r^', ms=15)
    if 'coj' in input_files[k]:
		plt.plot([mjd_obs[k] - mjd0], [np.power(np.sqrt(det_cov[k]),1./4.)], 'go', ms=15)
    #plt.plot(mjd_obs - mjd0, np.log10(np.array(det_cov)), 'ro')
    plot([mjd_obs[k] - mjd0,mjd_obs[k] - mjd0], [0.01,0.5], 'k--')
    plt.xlim(0.,2.3)
    plt.ylim(0.01,0.5)
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
    ylabel('$(\sqrt{|\Sigma|})^{1/4}$, mag')
    xlabel('Days Since MJD %1.0f'%mjd0)
    yticks([0.01, 0.03, 0.10, 0.3])

    ax=subplot(2,2,3)
    errorbar(mjd_obs - mjd0, mag1-1.0, yerr=(me1Lower, me1Upper), fmt=',', color=[0.5,0.5,0.5])
    errorbar([mjd_obs[k] - mjd0], [mag1[k]-1.0], yerr=([me1Lower[k]], [me1Upper[k]]), fmt=',', color=[0.,0.,0.],elinewidth=10)  #, label='M1 - 1.0'

    errorbar(mjd_obs - mjd0, mag3-0.5, yerr=(me3Lower, me3Upper), fmt=',', color=[0.2,0.2,1.])
    errorbar([mjd_obs[k] - mjd0], [mag3[k]-0.5], yerr=([me3Lower[k]], [me3Upper[k]]), fmt=',', color=[0.,0.,0.8],elinewidth=10) #, label='M2 - 0.5'

    errorbar(mjd_obs - mjd0, mag2, yerr=(me2Lower, me2Upper), fmt=',', color=[1.,0.2,0.2])
    errorbar([mjd_obs[k] - mjd0], [mag2[k]], yerr=([me2Lower[k]], [me2Upper[k]]), fmt=',', color=[0.7,0.,0.],elinewidth=10)     #, label='S1 + 0.0'

    errorbar(mjd_obs - mjd0, mag4+0.5, yerr=(me4Lower, me4Upper), fmt=',', color=[0.,1.,0.])
    errorbar([mjd_obs[k] - mjd0], [mag4[k]+0.5], yerr=([me4Lower[k]], [me4Upper[k]]), fmt=',', color='g',elinewidth=10) #, label='S2 + 0.5'

    plt.errorbar([-1.], [-1.], yerr=[1.], fmt='k,', elinewidth=10, label='M1 - 1.0')
    plt.errorbar([-1.], [-1.], yerr=[1.], fmt='b,',elinewidth=10,  label='M2 - 0.5')
    plt.errorbar([-1.], [-1.], yerr=[1.], fmt='r,',elinewidth=10,  label='S1 + 0.0')
    plt.errorbar([-1.], [-1.], yerr=[1.], fmt='g,',elinewidth=10,  label='S2 + 0.5')


    plot([mjd_obs[k] - mjd0,mjd_obs[k] - mjd0], [20.3,16.], 'k--')
    plt.xlim(0.,2.3)
    plt.ylim(20.3,16.)
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
    xlabel('Days Since MJD %1.0f'%mjd0)
    ylabel('Image Magnitudes')
    legend(loc=2, numpoints=1)
    suptitle('HE0435-1223 \n LCOGT Observations', fontsize=48)
    #suptitle(fnm.split('/')[-1])


    img_index = str(outFileTags[k]).split('_')[-1].split('.')[0]
    plt.savefig('slow_movie_image_%s.png'%img_index, dpi=50)
    plt.savefig('slow_movie_image_%s.pdf'%img_index)
    #plt.savefig('slow_movie_%s.png'%outFileTags[k], dpi=50)
    FM.hdulist.close()
    results.close()
    del FM
    print 'loop iteration finished: %1.2f s'%(clock()-start_time)




show()


import matplotlib
matplotlib.use('Agg')
import numpy as np
from pylab import *
from astropy.io import ascii
import glob
import AnalysisLibrary as AL
from time import clock

rcParams['font.size']=24
rcParams['legend.fontsize']=24
rcParams['figure.facecolor']='white'
#fnames = glob.glob('/nisushome/romerowo/lcolens_20150605/python/arw/npzfiles/image_*_results.npz') # no priors of seeing imposed
#fnames = glob.glob('/disk4/romerowo/lcolens_outputs/20150916/npzfiles/image_*_results.npz') # priors of seeing imposed
#fnames = glob.glob('/disk4/romerowo/lcolens_outputs/20150917/npzfiles/image_*_results.npz') # added pixelization and robust ZP (botched)
fnames = glob.glob('/disk4/romerowo/lcolens_outputs/20150918/npzfiles/image_*_results.npz') # robust ZP, no subsample-integrated pixelization
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
#rand_vals = arange(np.random.randint(0,len(fnames))-5,5)
#rand_vals = np.random.randint(0,len(fnames),5)
#rand_vals = sort(rand_vals)
#for k in [30, 80,90,100, 150, 210]:
#for k in rand_vals:
for k in range(0,len(fnames)):
	start_time = clock()
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
        #if results['chiSq'] > 2. or results['maxChi'] > 5.:
        if (1==1):
		if(k%20==0): print 'outFnm\t\tmjd_obs\t\tchiSq\tmaxChi\tinputFile\t\t\t\tm1\tZP_mean\tfilter'
		print '%s\t%1.3f\t%1.2f\t%1.2f\t%s\t%1.2f\t%1.2f\t%s'%(results['outFileTag'], results['mjd_obs'], results['chiSq'], results['maxChi'], str(results['inputFile']).split('/')[-1], results['m1'], results['ZP_mean'], results['filter'])
        input_files.append(str(results['inputFile']).split('/')[-1])
	if k==0:
		print '\t',results.files
	fnm = str(results['inputFile'])
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
        print 'esimating magnitudes: %1.2f s'%(clock()-start_time)
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

    	figure(figsize=(19.,24))

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
        cov = zeros((4,4))
        cov[0][0] = np.mean((mag1_array - mu1)**2)
        cov[1][1] = np.mean((mag2_array - mu2)**2)
        cov[2][2] = np.mean((mag3_array - mu3)**2)
        cov[3][3] = np.mean((mag4_array - mu4)**2)
        cov[0][1] = np.mean((mag1_array - mu1)*(mag2_array - mu2))
        cov[0][2] = np.mean((mag1_array - mu1)*(mag3_array - mu3))
        cov[0][3] = np.mean((mag1_array - mu1)*(mag4_array - mu4))
        cov[1][0] = cov[0][1]
        cov[1][2] = np.mean((mag2_array - mu2)*(mag3_array - mu3))
        cov[1][3] = np.mean((mag2_array - mu2)*(mag4_array - mu4))
        cov[2][0] = cov[0][2]
        cov[2][1] = cov[1][2]
        cov[2][3] = np.mean((mag3_array - mu3)*(mag4_array - mu4))
        cov[3][0] = cov[0][3]
        cov[3][1] = cov[1][3]
        cov[3][2] = cov[2][3]
        cov_det = np.linalg.det(cov)
        print '\t','***** COV DET **** %1.2e'%cov_det
	icov = np.infty
        icov_det = np.infty
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
    	ax1=plt.subplot(322)
    	ax2=plt.subplot(324)
    	ax3=plt.subplot(326)
	ZP_flx = AL.magnitude2flux(results['ZP_mean'])
        theta = [np.mean(x0Ch),np.mean(y0Ch),np.mean(amp1Ch),np.mean(amp2Ch),np.mean(amp3Ch),np.mean(amp4Ch), np.mean(alphaCh), np.mean(betaCh), np.mean(nbkgCh)]
        print 'plotting figures: %1.2f s'%(clock()-start_time)

	FM.plot_image_movie(ra_qsr, dec_qsr, ax1, ax2, ax3, ZP_flx, theta, Npx=31)
        plt.subplots_adjust(left=0.09)

        #print '\tPIXSCALE',pxscl
        ax=subplot(6,2,1)
        print '\t*********', results['inputFile']
        for nn in range(0,len(mjd_obs)):
            if 'lsc' in input_files[nn]:
	       plt.plot([mjd_obs[nn] - mjd0], [airmass[nn]], 'bs')
            if 'coj' in input_files[nn]:
	       plt.plot([mjd_obs[nn] - mjd0], [airmass[nn]], 'go')
            if 'cpt' in input_files[nn]:
	       plt.plot([mjd_obs[nn] - mjd0], [airmass[nn]], 'r^')
        plt.plot([-1.], [-1.], 'bs', ms=8, label='Chile')
        plt.plot([-1.], [-1.], 'go', ms=8, label='Australia')
        plt.plot([-1.], [-1.], 'r^', ms=8, label='S. Africa')
        plt.legend(loc=(0.17,0.), numpoints=1, fontsize=24, frameon=False, borderaxespad=-0.05, handletextpad=-0.5, columnspacing=-0.5, labelspacing=0., borderpad=-0.05)
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
	#plt.plot(mjd_obs - mjd0, seeing_fwhm, 'ro')
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
	#plt.plot(mjd_obs - mjd0, np.log10(np.array(det_cov)), 'ro')
        plt.xlim(0.,2.3)
        plt.ylim(0.01,0.5)
        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
	ylabel('$(\sqrt{|\Sigma|})^{1/4}$, mag')
        xlabel('Days Since MJD %1.0f'%mjd0)
        yticks([0.01, 0.03, 0.10, 0.3])

	'''
        ax=subplot(6,2,5)
	plt.plot(mjd_obs - mjd0, chiSq, 'ro')
        plt.xlim(0.,2.3)
        plt.ylim(0.,10.0)
        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
	ylabel('Chi Squared')
        xlabel('Days Since MJD %1.0f'%mjd0)
	'''

        ax=subplot(2,2,3)
	errorbar(mjd_obs - mjd0, mag1-1.0, yerr=(me1Lower, me1Upper), fmt=',', color=[0.5,0.5,0.5])
	errorbar(mjd_obs[-1:] - mjd0, mag1[-1:]-1.0, yerr=(me1Lower[-1:], me1Upper[-1:]), fmt=',', color='k',elinewidth=3, label='M1 - 1.0')
 
	errorbar(mjd_obs - mjd0, mag2-0.5, yerr=(me2Lower, me2Upper), fmt=',', color=[0.2,0.2,0.8])
	errorbar(mjd_obs[-1:] - mjd0, mag2[-1:]-0.5, yerr=(me2Lower[-1:], me2Upper[-1:]), fmt=',', color='b',elinewidth=3, label='M2 - 0.5')
 
	errorbar(mjd_obs - mjd0, mag3, yerr=(me3Lower, me3Upper), fmt=',', color=[0.8,0.2,0.2])
	errorbar(mjd_obs[-1:] - mjd0, mag3[-1:], yerr=(me3Lower[-1:], me3Upper[-1:]), fmt=',', color='r',elinewidth=3, label='S1 + 0.0')
 
	errorbar(mjd_obs - mjd0, mag4+0.5, yerr=(me4Lower, me4Upper), fmt=',', color=[0.2,0.8,0.2])
	errorbar(mjd_obs[-1:] - mjd0, mag4[-1:]+0.5, yerr=(me4Lower[-1:], me4Upper[-1:]), fmt=',', color='g',elinewidth=3, label='S2 + 0.5')
        plt.xlim(0.,2.3)
        plt.ylim(20.3,16.)
        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
        xlabel('Days Since MJD %1.0f'%mjd0)
	ylabel('Image Magnitudes')
	legend(loc=2, numpoints=1)
	suptitle('HE0435-1223 \n LCOGT Observations', fontsize=48)
	#suptitle(fnm.split('/')[-1])

	'''
        ax=subplot(2,2,3)
	errorbar(mjd_obs - mjd0, flux1, yerr=(flux1LowerError, flux1UpperError), fmt=',', color='k')
	errorbar(mjd_obs[-1:] - mjd0, flux1[-1:], yerr=(flux1LowerError[-1:], flux1UpperError[-1:]), fmt=',', color='red',elinewidth=3)
        plt.xlim(0.,2.3)
        plt.ylim(0.5e-8,5.5e-8)
        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
	ylabel('Flux Image 1')

	errorbar(mjd_obs - mjd0, flux2, yerr=(flux2LowerError, flux2UpperError), fmt=',', color='k')
	errorbar(mjd_obs[-1:] - mjd0, flux2[-1:], yerr=(flux2LowerError[-1:], flux2UpperError[-1:]), fmt=',', color='red',elinewidth=3)
        plt.xlim(0.,2.3)
        plt.ylim(0.5e-8,5.5e-8)
        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
	ylabel('Flux Image 2')

	errorbar(mjd_obs - mjd0, flux3, yerr=(flux2LowerError, flux2UpperError), fmt=',', color='k')
	errorbar(mjd_obs[-1:] - mjd0, flux3[-1:], yerr=(flux3LowerError[-1:], flux3UpperError[-1:]), fmt=',', color='red',elinewidth=3)
        plt.xlim(0.,2.3)
        plt.ylim(0.5e-8,5.5e-8)
        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
	ylabel('Flux Image 3')

	errorbar(mjd_obs - mjd0, flux4, yerr=(flux2LowerError, flux2UpperError), fmt=',', color='k')
	errorbar(mjd_obs[-1:] - mjd0, flux4[-1:], yerr=(flux1LowerError[-1:], flux4UpperError[-1:]), fmt=',', color='red',elinewidth=3)
        plt.xlim(0.,2.3)
        plt.ylim(0.5e-8,5.5e-8)
        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
        xlabel('Days Since MJD %1.0f'%mjd0)
	ylabel('Flux Image 4')
	'''

    	plt.savefig('movie_%s.png'%results['outFileTag'], dpi=50)
	FM.hdulist.close()
	results.close()
        print 'loop iteration finished: %1.2f s'%(clock()-start_time)


exit()

#hist(fwhmCh)

mjd_obs = array(mjd_obs)
m1 = array(m1)
m2 = array(m2)
m3 = array(m3)
m4 = array(m4)
me1 = array(me1)
me2 = array(me2)
me3 = array(me3)
me4 = array(me4)
err_array = np.concatenate([me1, me2, me3, me4])
maxChi=array(maxChi)
chiSq=array(chiSq)

ZP_mean = np.array(ZP_mean)
#ZP_std = np.array(ZP_std)
ZP_wrms = np.array(ZP_wrms)

amp1Emcee = np.array(amp1Emcee)
amp1EmceeUpperError = np.array(amp1EmceeUpperError)
amp1EmceeLowerError = np.array(amp1EmceeLowerError)
amp2Emcee = np.array(amp2Emcee)
amp2EmceeUpperError = np.array(amp2EmceeUpperError)
amp2EmceeLowerError = np.array(amp2EmceeLowerError)
amp3Emcee = np.array(amp3Emcee)
amp3EmceeUpperError = np.array(amp3EmceeUpperError)
amp3EmceeLowerError = np.array(amp3EmceeLowerError)
amp4Emcee = np.array(amp4Emcee)
amp4EmceeUpperError = np.array(amp4EmceeUpperError)
amp4EmceeLowerError = np.array(amp4EmceeLowerError)
nbkgEmcee = np.array(nbkgEmcee)
nbkgEmceeUpperError = np.array(nbkgEmceeUpperError)
nbkgEmceeLowerError = np.array(nbkgEmceeLowerError)

LC1Emcee = np.array(LC1Emcee)
LC1EmceeUpperError = np.array(LC1EmceeUpperError)
LC1EmceeLowerError = np.array(LC1EmceeLowerError)
LC2Emcee = np.array(LC2Emcee)
LC2EmceeUpperError = np.array(LC2EmceeUpperError)
LC2EmceeLowerError = np.array(LC2EmceeLowerError)
LC3Emcee = np.array(LC3Emcee)
LC3EmceeUpperError = np.array(LC3EmceeUpperError)
LC3EmceeLowerError = np.array(LC3EmceeLowerError)
LC4Emcee = np.array(LC4Emcee)
LC4EmceeUpperError = np.array(LC4EmceeUpperError)
LC4EmceeLowerError = np.array(LC4EmceeLowerError)
#figure()
#plot(maxChi, '.')
err_cut = 0.15
cut12   = np.logical_and(me2<err_cut, me1<err_cut)
cut123  = np.logical_and(me3<err_cut, cut12)
cut1234 = np.logical_and(me4<err_cut, cut123)
cut = np.logical_and(maxChi<5.,cut1234)

print len(mjd_obs)
print len(mjd_obs[cut])
ascii.write([mjd_obs[cut], m1[cut], me1[cut], m2[cut], me2[cut], m3[cut], me3[cut], m4[cut], me4[cut]], 'he0435-1223_lcogt_magnitudes.dat', names=['mjd', 'mag_A', 'magerr_A', 'mag_B', 'magerr_B', 'mag_C', 'magerr_C', 'mag_D', 'magerr_D'] )

mjd_obs = array(mjd_obs)
mjd0 = np.floor(np.min(mjd_obs))

mx = max([max(m1), max(m2), max(m3), max(m4)])+0.2
mn = min([min(m1), min(m2), min(m3), min(m4)])-0.2

figure()
for k in range(0,len(site)):

	if(site[k]=='lsc'):
		p1, = plot(mjd_obs[k] - mjd0, readnoise[k], 'bo', label='lsc')
		#errorbar(mjd_obs[k] - mjd0, readnoise[k], yerr=readnoise_error[k], fmt='b,', alpha=0.5)
	if(site[k]=='cpt'):
		p2, = plot(mjd_obs[k] - mjd0, readnoise[k], 'go', label='cpt')
		#errorbar(mjd_obs[k] - mjd0, readnoise[k], yerr=readnoise_error[k], fmt='g,', alpha=0.5)
	if(site[k]=='coj'):
		p3, = plot(mjd_obs[k] - mjd0, readnoise[k], 'ro', label='coj')
		#errorbar(mjd_obs[k] - mjd0, readnoise[k], yerr=readnoise_error[k], fmt='r,', alpha=0.5)
try:
	legend((p1,p2,p3), ('lsc', 'cpt', 'coj'), loc=4, numpoints=1)
except:
	print 'nevermind'
xlabel('Days Since mjd %1.0f'%(mjd0))
ylabel('Estimated Read Noise, counts')
title('Read Noise by Observation')
ylim(0.,35.)
grid(True)

figure()
fwhm =  array(alpha)*2.*np.sqrt(2.**(1/array(beta))-1.)
#plot(mjd_obs[k] - mjd0, fwhm[k], 'ro', label='coj')
for k in range(0,len(site)):
	if(site[k]=='lsc'):
		p1, = plot(mjd_obs[k] - mjd0, fwhm[k], 'bo', label='lsc')
	if(site[k]=='cpt'):
		p2, = plot(mjd_obs[k] - mjd0, fwhm[k], 'go', label='cpt')
	if(site[k]=='coj'):
		p3, = plot(mjd_obs[k] - mjd0, fwhm[k], 'ro', label='coj')
try:
	legend((p1,p2,p3), ('lsc', 'cpt', 'coj'), loc=0, numpoints=1)
except:
	print 'nevermind'
xlabel('Days Since mjd %1.0f'%(mjd0))
ylabel('FWHM, pixels')
title('Moffat FWHM')
grid(True)
ylim(0.,5.5)

figure()
for k in range(0,len(site)):
	if(site[k]=='lsc'):
		#print fwhmEmceeLowerError[k],fwhmEmceeUpperError[k]
		p1, = plot(mjd_obs[k] - mjd0, fwhmEmcee[k], 'b*', label='lsc')
		errorbar(mjd_obs[k] - mjd0, fwhmEmcee[k], yerr=[[fwhmEmceeLowerError[k]],[fwhmEmceeUpperError[k]]], fmt='b,')
	if(site[k]=='cpt'):
		p2, = plot(mjd_obs[k] - mjd0, fwhmEmcee[k], 'g*', label='cpt')
		errorbar(mjd_obs[k] - mjd0, fwhmEmcee[k], yerr=[[fwhmEmceeLowerError[k]],[fwhmEmceeUpperError[k]]], fmt='g,')
	if(site[k]=='coj'):
		p3, = plot(mjd_obs[k] - mjd0, fwhmEmcee[k], 'r*', label='coj')
		errorbar(mjd_obs[k] - mjd0, fwhmEmcee[k], yerr=[[fwhmEmceeLowerError[k]],[fwhmEmceeUpperError[k]]], fmt='r,')
try:
	legend((p1,p2,p3), ('lsc', 'cpt', 'coj'), loc=0, numpoints=1)
except:
	print 'nevermind'
xlabel('Days Since mjd %1.0f'%(mjd0))
ylabel('FWHM, pixels')
title('Moffat FWHM')
grid(True)
ylim(0.,5.5)


#plot(np.ma.masked_where(site=='lsc',mjd_obs) - mjd0, np.ma.masked_where(site=='lsc',readnoise), 'bo')
#plot(np.ma.masked_where(site=='coj',mjd_obs) - mjd0, np.ma.masked_where(site=='coj',readnoise), 'ro')
#plot(np.ma.masked_where(site=='cpt',mjd_obs) - mjd0, np.ma.masked_where(site=='cpt',readnoise), 'go')

figure(figsize=(8,12))
'''
for k in range(0,len(mjd_obs)):

  col = 'r'
  if(filter[k]=='gp'):
  	col='g'
  subplot(311)
  errorbar(mjd_obs[k] - mjd0, ZP_mean[k], yerr=ZP_wrms[k],fmt='o', color=col)
  grid(True)
  #subplot(412)
  #plot(mjd_obs[k] - mjd0, ZP_rms[k], 'o', color=col)
  #grid(True)
  subplot(312)
  plot(mjd_obs[k] - mjd0, ZP_wrms[k], 'o', color=col)
  grid(True)
  subplot(313)
  plot(mjd_obs[k] - mjd0, ZP_std[k], 'o', color=col)  
'''
subplot(311)
gca().invert_yaxis()
ylabel('ZP Weighted Mean, mag.')
#subplot(412)
#ylabel('ZP RMS, mag.')
subplot(312)
ylabel('ZP wRMS, mag.')
subplot(313)
grid(True)
xlabel('Days Since mjd %1.0f'%(mjd0))
ylabel('Uncertainty on the\nZP Mean, mag.')
subplots_adjust(left=0.2)


figure(figsize=(8,12))
ax1 = subplot(411)
errorbar(mjd_obs[cut] - mjd0, amp1Emcee[cut], yerr=(amp1EmceeLowerError[cut], amp1EmceeUpperError[cut]), fmt=',', color='k')
grid(True, which="both")
ylabel('image 1\nMoffat Amplitude')
ax1 = subplot(412)
errorbar(mjd_obs[cut] - mjd0, amp2Emcee[cut], yerr=(amp2EmceeLowerError[cut], amp2EmceeUpperError[cut]), fmt=',', color='k')
grid(True, which="both")
ylabel('image 2\nMoffat Amplitude')
ax1 = subplot(413)
errorbar(mjd_obs[cut] - mjd0, amp3Emcee[cut], yerr=(amp3EmceeLowerError[cut], amp3EmceeUpperError[cut]), fmt=',', color='k')
grid(True, which="both")
ylabel('image 3\nMoffat Amplitude')
ax1 = subplot(414)
errorbar(mjd_obs[cut] - mjd0, amp4Emcee[cut], yerr=(amp4EmceeLowerError[cut], amp4EmceeUpperError[cut]), fmt=',', color='k')
grid(True, which="both")
ylabel('image 4\nMoffat Amplitude')
xlabel('Days Since mjd %1.0f'%(mjd0))

figure(figsize=(8,12))
ax1 = subplot(411)
errorbar(mjd_obs[cut] - mjd0, LC1Emcee[cut], yerr=(LC1EmceeLowerError[cut], LC1EmceeUpperError[cut]), fmt=',', color='k')
grid(True, which="both")
ylabel('image 1\nInstrument Flux')
ax1 = subplot(412)
errorbar(mjd_obs[cut] - mjd0, LC2Emcee[cut], yerr=(LC2EmceeLowerError[cut], LC2EmceeUpperError[cut]), fmt=',', color='k')
grid(True, which="both")
ylabel('image 2\nInstrument Flux')
ax1 = subplot(413)
errorbar(mjd_obs[cut] - mjd0, LC3Emcee[cut], yerr=(LC3EmceeLowerError[cut], LC3EmceeUpperError[cut]), fmt=',', color='k')
grid(True, which="both")
ylabel('image 3\nInstrument Flux')
ax1 = subplot(414)
errorbar(mjd_obs[cut] - mjd0, LC4Emcee[cut], yerr=(LC4EmceeLowerError[cut], LC4EmceeUpperError[cut]), fmt=',', color='k')
grid(True, which="both")
ylabel('image 4\nInstrument Flux')
xlabel('Days Since mjd %1.0f'%(mjd0))
subplots_adjust(left=0.2)

ZP_flux       = AL.magnitude2flux(ZP_mean)
#ZP_flux_std   = AL.magErr2fluxErr(ZP_mean, ZP_std)
ZP_flux_wrms  = AL.magErr2fluxErr(ZP_mean, ZP_wrms)
flux1 = np.array(ZP_flux*LC1Emcee)
flux1UpperError=np.array(flux1*np.sqrt(LC1EmceeUpperError**2/LC1Emcee**2 + ZP_flux_wrms**2/ZP_flux**2))
flux1LowerError=np.array(flux1*np.sqrt(LC1EmceeLowerError**2/LC1Emcee**2 + ZP_flux_wrms**2/ZP_flux**2))
flux2 = np.array(ZP_flux*LC2Emcee)
flux2UpperError=np.array(flux2*np.sqrt(LC2EmceeUpperError**2/LC2Emcee**2 + ZP_flux_wrms**2/ZP_flux**2))
flux2LowerError=np.array(flux2*np.sqrt(LC2EmceeLowerError**2/LC2Emcee**2 + ZP_flux_wrms**2/ZP_flux**2))
flux3 = np.array(ZP_flux*LC3Emcee)
flux3UpperError=np.array(flux3*np.sqrt(LC3EmceeUpperError**2/LC3Emcee**2 + ZP_flux_wrms**2/ZP_flux**2))
flux3LowerError=np.array(flux3*np.sqrt(LC3EmceeLowerError**2/LC3Emcee**2 + ZP_flux_wrms**2/ZP_flux**2))
flux4 = np.array(ZP_flux*LC4Emcee)
flux4UpperError=np.array(flux4*np.sqrt(LC4EmceeUpperError**2/LC4Emcee**2 + ZP_flux_wrms**2/ZP_flux**2))
flux4LowerError=np.array(flux4*np.sqrt(LC4EmceeLowerError**2/LC4Emcee**2 + ZP_flux_wrms**2/ZP_flux**2))

figure(figsize=(8,12))
ax1 = subplot(411)
errorbar(mjd_obs[cut] - mjd0, flux1[cut], yerr=(flux1LowerError[cut], flux1UpperError[cut]), fmt=',', color='k')
grid(True, which="both")
ylabel('image 1 Flux')
ax1 = subplot(412)
errorbar(mjd_obs[cut] - mjd0, flux2[cut], yerr=(flux2LowerError[cut], flux2UpperError[cut]), fmt=',', color='k')
grid(True, which="both")
ylabel('image 2  Flux')
ax1 = subplot(413)
errorbar(mjd_obs[cut] - mjd0, flux3[cut], yerr=(flux3LowerError[cut], flux3UpperError[cut]), fmt=',', color='k')
grid(True, which="both")
ylabel('image 3  Flux')
ax1 = subplot(414)
errorbar(mjd_obs[cut] - mjd0, flux4[cut], yerr=(flux4LowerError[cut], flux4UpperError[cut]), fmt=',', color='k')
grid(True, which="both")
ylabel('image 4 Flux')
xlabel('Days Since mjd %1.0f'%(mjd0))
subplots_adjust(left=0.2)

m1 = AL.flux2magnitude(flux1)
me1Upper = AL.fluxErr2magErr(flux1, flux1UpperError)
me1Lower = AL.fluxErr2magErr(flux1, flux1LowerError)
m2 = AL.flux2magnitude(flux2)
me2Upper = AL.fluxErr2magErr(flux2, flux2UpperError)
me2Lower = AL.fluxErr2magErr(flux2, flux2LowerError)
m3 = AL.flux2magnitude(flux3)
me3Upper = AL.fluxErr2magErr(flux3, flux3UpperError)
me3Lower = AL.fluxErr2magErr(flux3, flux3LowerError)
m4 = AL.flux2magnitude(flux4)
me4Upper = AL.fluxErr2magErr(flux4, flux4UpperError)
me4Lower = AL.fluxErr2magErr(flux4, flux4LowerError)

figure(figsize=(8,12))
ax1 = subplot(411)
errorbar(mjd_obs[cut] - mjd0, m1[cut], yerr=(me1Lower[cut], me1Upper[cut]), fmt=',', color='k')
gca().invert_yaxis()
grid(True, which="both")
ylabel('image 1 mag.')
ax1 = subplot(412)
errorbar(mjd_obs[cut] - mjd0, m2[cut], yerr=(me2Lower[cut], me2Upper[cut]), fmt=',', color='k')
gca().invert_yaxis()
grid(True, which="both")
ylabel('image 2  mag.')
ax1 = subplot(413)
errorbar(mjd_obs[cut] - mjd0, m3[cut], yerr=(me3Lower[cut], me3Upper[cut]), fmt=',', color='k')
gca().invert_yaxis()
grid(True, which="both")
ylabel('image 3  mag.')
ax1 = subplot(414)
errorbar(mjd_obs[cut] - mjd0, m4[cut], yerr=(me4Lower[cut], me4Upper[cut]), fmt=',', color='k')
gca().invert_yaxis()
grid(True, which="both")
ylabel('image 4 mag.')
xlabel('Days Since mjd %1.0f'%(mjd0))
subplots_adjust(left=0.2)

figure(figsize=(8,12))
ax1 = subplot(611)
semilogy(mjd_obs[cut] - mjd0, ZP_flux_std[cut]/ZP_flux[cut], '.', color='r')
ylabel('zero point\nerror on the mean')
grid(True, which="both")
ylim(1.e-3, 1.0)
title('Factional Errors on Flux')
ax1 = subplot(612)
semilogy(mjd_obs[cut] - mjd0, ZP_flux_wrms[cut]/ZP_flux[cut], '.', color='r')
ylabel('zero point\nwrms')
grid(True, which="both")
ylim(1.e-3, 1.0)
ax1 = subplot(613)
semilogy(mjd_obs[cut] - mjd0, 0.5*(flux1LowerError[cut] + flux1UpperError[cut])/flux1[cut], '.', color='r')
grid(True, which="both")
ylabel('image 1')
ylim(1.e-3, 1.0)
ax1 = subplot(614)
semilogy(mjd_obs[cut] - mjd0, 0.5*(flux2LowerError[cut] + flux2UpperError[cut])/flux2[cut], '.', color='r')
grid(True, which="both")
ylabel('image 2')
ylim(1.e-3, 1.0)
ax1 = subplot(615)
semilogy(mjd_obs[cut] - mjd0, 0.5*(flux3LowerError[cut] + flux3UpperError[cut])/flux3[cut], '.', color='r')
grid(True, which="both")
ylabel('image 3')
ylim(1.e-3, 1.0)
ax1 = subplot(616)
semilogy(mjd_obs[cut] - mjd0, 0.5*(flux4LowerError[cut] + flux4UpperError[cut])/flux4[cut], '.', color='r')
grid(True, which="both")
ylabel('image 4')
xlabel('Days Since mjd %1.0f'%(mjd0))
ylim(1.e-3, 1.0)
subplots_adjust(left=0.2)


figure(figsize=(8,8))
semilogy(mjd_obs[cut] - mjd0, ZP_flux_std[cut]/ZP_flux[cut], 'o', label='zero point\nerror on the mean')
title('Factional Errors on Flux')
semilogy(mjd_obs[cut] - mjd0, ZP_flux_wrms[cut]/ZP_flux[cut], 's', label='zero point\nwrms')
semilogy(mjd_obs[cut] - mjd0, 0.5*(flux1LowerError[cut] + flux1UpperError[cut])/flux1[cut], '<', label='image 1')
semilogy(mjd_obs[cut] - mjd0, 0.5*(flux2LowerError[cut] + flux2UpperError[cut])/flux2[cut], '>', label='image 2')
semilogy(mjd_obs[cut] - mjd0, 0.5*(flux3LowerError[cut] + flux3UpperError[cut])/flux3[cut], '^', label='image 3')
semilogy(mjd_obs[cut] - mjd0, 0.5*(flux4LowerError[cut] + flux4UpperError[cut])/flux4[cut], 'v', label='image 4')
grid(True, which="both")
xlabel('Days Since mjd %1.0f'%(mjd0))
ylim(1.e-3, 1.0)
legend(loc=1, numpoints=1)
subplots_adjust(left=0.2)

show()

figure(figsize=(8,18))
ax1 = subplot(511)
#errorbar(mjd_obs - mjd0, m1, yerr=me1, fmt='o', color='k')
errorbar(mjd_obs[cut] - mjd0, m1[cut], yerr=me1[cut], fmt='o', color='k')
gca().invert_yaxis()
grid(True, which="both")
ylabel('image 1 magnitude')
ylim(30.,0.)
#ylim(mx,mn)

subplot(512, sharex=ax1, sharey=ax1)
#errorbar(mjd_obs - mjd0, m2, yerr=me2, fmt='ko', color='b')
errorbar(mjd_obs[cut] - mjd0, m2[cut], yerr=me2[cut], fmt='ko', color='b')
gca().invert_yaxis()
ylim(30.,0.)
#ylim(mx,mn)
grid(True, which="both")
ylabel('image 2 magnitude')

subplot(513, sharex=ax1, sharey=ax1)
#errorbar(mjd_obs - mjd0, m3, yerr=me3, fmt='ko', color='r')
errorbar(mjd_obs[cut] - mjd0, m3[cut], yerr=me3[cut], fmt='ko', color='r')
gca().invert_yaxis()
ylim(30.,0.)
#ylim(mx,mn)
grid(True, which="both")
ylabel('image 3 magnitude')

subplot(514, sharex=ax1, sharey=ax1)
#errorbar(mjd_obs - mjd0, m4, yerr=me4, fmt='ko', color='g')
errorbar(mjd_obs[cut]- mjd0, m4[cut], yerr=me4[cut], fmt='ko', color='g')
gca().invert_yaxis()
ylim(30.,0.)
#ylim(mx,mn)
grid(True, which="both")
ylabel('image 4 magnitude')

ax=subplot(515, sharex=ax1)
ax.set_yscale('log')
#plot(mjd_obs - mjd0, chiSq, 'bo', label='Chi Square')
#plot(mjd_obs - mjd0, maxChi, 'ro', label='Max. Chi Value')
plot(mjd_obs[cut] - mjd0, chiSq[cut], 'bo', label='Chi Square')
plot(mjd_obs[cut] - mjd0, maxChi[cut], 'ro', label='Max. Chi Value')
legend(loc=0, title='Quad. Fit Image', numpoints=1)
xlabel('Days Since mjd %1.0f'%(mjd0))
grid(True, which="both")
suptitle('HE0435-1223 LCOGT Light Curves')
subplots_adjust(hspace=0.3, bottom=0.05, top=0.95)
ylim(0.1,1000.)
#subplot(615)

#subplot(616)


show()

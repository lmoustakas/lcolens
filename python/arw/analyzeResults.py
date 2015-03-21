import numpy as np
from pylab import *
from astropy.io import ascii
import glob

rcParams['font.size']=14
rcParams['legend.fontsize']=14
rcParams['figure.facecolor']='white'
fnames = glob.glob('image_*_results.npz')
fnames = sorted(fnames)
for f in fnames:
	print f

#results = np.load('out_results.npz')
#print npzfile.files
#print npzfile['m1'], npzfile['me1']


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
site=[]
ZP_mean=[]
ZP_rms=[]
filter=[]
alpha=[]
beta=[]
for k in range(0,len(fnames)):
	results = np.load(fnames[k])
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
	site.append(str(results['inputFile']).split('/')[-1][0:3])
	ZP_mean.append(results['ZP_mean'])
	ZP_rms.append(results['ZP_rms'])
	filter.append(results['filter'])
	alpha.append(results['APASS_alpha'])
	beta.append(results['APASS_beta'])
        #if results['chiSq'] > 2. or results['maxChi'] > 5.:
        if (1==1):
		if(k%20==0): print 'outFnm\t\tmjd_obs\t\tchiSq\tmaxChi\tinputFile\t\t\t\tm1\tZP_mean\tfilter'
		print '%s\t%1.3f\t%1.2f\t%1.2f\t%s\t%1.2f\t%1.2f\t%s'%(results['outFileTag'], results['mjd_obs'], results['chiSq'], results['maxChi'], str(results['inputFile']).split('/')[-1], results['m1'], results['ZP_mean'], results['filter'])
	if k==0:
		print results.files

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
#figure()
#plot(maxChi, '.')
#show()
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
	if(site[k]=='cpt'):
		p2, = plot(mjd_obs[k] - mjd0, readnoise[k], 'go', label='cpt')
	if(site[k]=='coj'):
		p3, = plot(mjd_obs[k] - mjd0, readnoise[k], 'ro', label='coj')
legend((p1,p2,p3), ('lsc', 'cpt', 'coj'), loc=4, numpoints=1)
xlabel('Days Since mjd %1.0f'%(mjd0))
ylabel('Estimated Read Noise, counts')
title('Read Noise by Observation')
grid(True)

figure()
fwhm =  array(alpha)*2.*np.sqrt(2.**(1/array(beta))-1.)
plot(mjd_obs[k] - mjd0, fwhm[k], 'ro', label='coj')
for k in range(0,len(site)):

	if(site[k]=='lsc'):
		p1, = plot(mjd_obs[k] - mjd0, fwhm[k], 'bo', label='lsc')
	if(site[k]=='cpt'):
		p2, = plot(mjd_obs[k] - mjd0, fwhm[k], 'go', label='cpt')
	if(site[k]=='coj'):
		p3, = plot(mjd_obs[k] - mjd0, fwhm[k], 'ro', label='coj')
legend((p1,p2,p3), ('lsc', 'cpt', 'coj'), loc=0, numpoints=1)
xlabel('Days Since mjd %1.0f'%(mjd0))
ylabel('FWHM, pixels')
title('Moffat FWHM')
grid(True)
#plot(np.ma.masked_where(site=='lsc',mjd_obs) - mjd0, np.ma.masked_where(site=='lsc',readnoise), 'bo')
#plot(np.ma.masked_where(site=='coj',mjd_obs) - mjd0, np.ma.masked_where(site=='coj',readnoise), 'ro')
#plot(np.ma.masked_where(site=='cpt',mjd_obs) - mjd0, np.ma.masked_where(site=='cpt',readnoise), 'go')
#show()

figure()
for k in range(0,len(mjd_obs)):
  if(filter[k]=='rp'):
    errorbar(mjd_obs[k] - mjd0, ZP_mean[k], yerr=ZP_rms[k], fmt='o', color='r')
  if(filter[k]=='gp'):
    errorbar(mjd_obs[k] - mjd0, ZP_mean[k], yerr=ZP_rms[k], fmt='o', color='g')
grid(True)
xlabel('Days Since mjd %1.0f'%(mjd0))
ylabel('Zero Point, magnitude')

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

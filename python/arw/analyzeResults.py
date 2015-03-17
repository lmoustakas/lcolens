import numpy as np
from pylab import *
import glob

rcParams['font.size']=14
rcParams['legend.fontsize']=12
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
        if results['chiSq'] > 2. or results['maxChi'] > 5.:
		print '%s\t%1.2f\t%1.2f\t%s'%(results['outFileTag'], results['chiSq'], results['maxChi'], str(results['inputFile']).split('/')[-1])
	if k==0:
		print results.files

mjd_obs = array(mjd_obs)
mjd0 = np.floor(np.min(mjd_obs))

mx = max([max(m1), max(m2), max(m3), max(m4)])+0.2
mn = min([min(m1), min(m2), min(m3), min(m4)])-0.2
figure(figsize=(8,18))
ax1 = subplot(511)
errorbar(mjd_obs - mjd0, m1, yerr=me1, fmt='o', color='k')
gca().invert_yaxis()
grid(True, which="both")
ylabel('image 1 magnitude')

ylim(mx,mn)
subplot(512, sharex=ax1)
errorbar(mjd_obs - mjd0, m2, yerr=me2, fmt='ko', color='b')
gca().invert_yaxis()
ylim(mx,mn)
grid(True, which="both")
ylabel('image 2 magnitude')

subplot(513, sharex=ax1)
errorbar(mjd_obs - mjd0, m3, yerr=me3, fmt='ko', color='r')
gca().invert_yaxis()
ylim(mx,mn)
grid(True, which="both")
ylabel('image 3 magnitude')

subplot(514, sharex=ax1)
errorbar(mjd_obs - mjd0, m4, yerr=me4, fmt='ko', color='g')
gca().invert_yaxis()
ylim(mx,mn)
grid(True, which="both")
ylabel('image 4 magnitude')

ax=subplot(515, sharex=ax1)
ax.set_yscale('log')
plot(mjd_obs - mjd0, chiSq, 'bo', label='Chi Square')
plot(mjd_obs - mjd0, maxChi, 'ro', label='Max. Chi Value')
legend(loc=0, title='Quad. Fit Image', numpoints=1)
xlabel('Days Since mjd %1.0f'%(mjd0))
grid(True, which="both")
suptitle('HE0435-1223 LCOGT Light Curves')
subplots_adjust(hspace=0.3, bottom=0.05, top=0.95)
#subplot(615)

#subplot(616)


show()

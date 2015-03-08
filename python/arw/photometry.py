from AnalysisLibrary import *
from pylab import *
import datetime

rcParams['figure.facecolor']='white'

# SUPER CRUDE PHOTOMETRY

# 1. Read APASS Catalogue
# 2. Loop through data
# 3. Estimage N_bkg from distribution of counts for full FITS image
# 4. Loop through APASS sources
#    a. get the peak value of the image to 
#    b. compare magnitude to APASS magnitude
#    c. derive correction
#    d. look at mean and rms of correction.
#    e. get quasar image
#    f. get peak values of the four images
#    g. correct the values to magnitudes
#    h. derive time-series with root N errors.

dirc = '/data2/romerowo/lcogt_data/he045-1223_wcs_corrected/'

APASS_table = ascii.read('../../data/HE0435_LCOGT/APASS_0438_list.csv')
#filename_table = ascii.read('time_ordered_file_list.dat')
filename_table = ascii.read('t_ord_image_stats.dat')
print APASS_table

#ra = 04:38:14.9000 
ra = (4.+38./60.+14.9/60./60.)/24*360 
#dec = -12:17:14.40
dec = -12. - (17./60 +14.4/60./60.)

APASS_rejects = [9, 21, 22] # 9 and 21 seem to be at the edge of the field of view. 22 seems close to the edge, sometimes we don't catch it.

count = 0
rando = np.random.randint(37,400)
#for fnm in filename_table['filename']:
for k in range(0,len(filename_table['filename'])):
  count += 1
  fnm = filename_table['filename'][k]
  if(filename_table['mjd'][count-1]<57008.):
	continue
  FM = FITSmanager(dirc+fnm)
  #FM.histogram_image_values(NBINS=5000, rng_max=5000.)
  N_bkg = float(filename_table['N_bkg'][k])
  sigma_read = float(filename_table['read_noise'][k])
  print 'N_bkg', N_bkg
  # CAN ESTIMATE THE READ NOISE BY TAKING THE WIDTH OF THIS DISTRIBUTION
  #FM.plot_image_values(NBINS=10000, rng_max=5000.)
  #show()
  filt = 'Sloan_r'
  if( FM.hdulist[0].header['FILTER'] == 'gp'):
    filt = 'Sloan_g'
  t1 = FM.hdulist[0].header['UTSTART']
  t2 = FM.hdulist[0].header['UTSTOP']
  #print t1+'000',t2
  t1 = datetime.datetime.strptime(t1,"%H:%M:%S.%f")
  t2 = datetime.datetime.strptime(t2,"%H:%M:%S.%f")
  t = (t2-t1).seconds + (t2-t1).microseconds/1.e6
  print '%d'%(count), fnm
  v1 = []
  v2 = []
  
  for k in range(0,len(APASS_table)):
    if(APASS_table[filt][k]!='NA' and k not in APASS_rejects): # entry 21 is on the edge of the field of view.
	obj = SourceImage(FM, APASS_table['radeg'][k], APASS_table['decdeg'][k], 30)
	# SUBTRACT BACKGROUND COUNTS
	obj.image -= N_bkg
	if(len(obj.image)==30):
  	  #imshow(obj.image)
	  print '%d.%d'%(count,k), fnm, len(obj.image), APASS_table['radeg'][k], APASS_table['decdeg'][k]
	  m_I = -2.5*log10(np.max(obj.image)/t)
	  print '\t',filt, APASS_table[filt][k], m_I
	  v1.append(float(APASS_table[filt][k]))
	  v2.append(m_I)



  figure(1)
  v1,v2 = zip(*sorted(zip(v1,v2)))
  ZP = array(v1)-array(v2)
  if(filt=='Sloan_r'):
	subplot(211)
	plot(v1,array(v1)-array(v2),'.-')
	#subplot(222)
	#hist(array(v1)-array(v2))
  if(filt=='Sloan_g'):
	subplot(212)
	plot(v1,array(v1)-array(v2),'.-')
	#subplot(224)
	#hist(array(v1)-array(v2))

  figure(2)
  ZP_mean = sum(ZP)/len(ZP)
  ZP_rms = sqrt(sum(ZP**2)/len(ZP) - ZP_mean**2)
  print '\t', ZP_mean, ZP_rms
  if(filt=='Sloan_r'):
	errorbar([filename_table['mjd'][count-1]-57008.],[ZP_mean],yerr=[ZP_rms],fmt='ro')
  if(filt=='Sloan_g'):
	errorbar([filename_table['mjd'][count-1]-57008.],[ZP_mean],yerr=[ZP_rms],fmt='go')

  obj = SourceImage(FM, ra, dec, 30)

  figure(10, figsize=(7,12))
  col = 'red'
  if(filt=='Sloan_g'): col='green'
  subplot(511)
  errorbar([filename_table['mjd'][count-1]-57008.],[N_bkg],yerr=[sqrt(N_bkg+sigma_read**2)],fmt='s', color=col, label='QSR peak')
  errorbar([filename_table['mjd'][count-1]-57008.],[np.max(obj.image)],yerr=[sqrt(np.max(obj.image)+sigma_read**2)],fmt='*', color=col, label='Background')
  ylabel('Raw Counts')
  #legend(loc=0, numpoints=1)

  subplot(512)
  errorbar([filename_table['mjd'][count-1]-57008.],[np.max(obj.image)-N_bkg],yerr=[sqrt(np.max(obj.image)+sigma_read**2)],fmt='o', color=col)
  ylabel('$N_{QSR}-N_{bkg}$,\nRaw Counts')

  obj.image -= N_bkg
  pk1 = np.max(obj.image)
  print '\t', pk1, -2.5*log10(pk1/t),  -2.5*log10(pk1/t) + ZP_mean 
  err1 = sqrt(pk1 + (sigma_read/t)**2)
  m_pk = -2.5*log10(pk1/t)
  m_err1 =  2.5/err1/log(10.)
  m_err1 = sqrt(m_err1**2)

  subplot(513)
  errorbar([filename_table['mjd'][count-1]-57008.],[m_pk],yerr=[m_err1],fmt='o', color=col)
  ylabel(r'$m_{I}=-2.5 log_{10}((N_{QSR}-N_{bkg})/t)$')
  gca().invert_yaxis()
  m_pk = -2.5*log10(pk1/t) + ZP_mean 
  m_err1 =  2.5/err1/log(10.)
  m_err1 = sqrt(m_err1**2 + ZP_rms**2 )

  subplot(514)
  errorbar([filename_table['mjd'][count-1]-57008.],[ZP_mean],yerr=[ZP_rms],fmt='o', color=col)

  ylabel(r'Zero Point')

  subplot(515)
  errorbar([filename_table['mjd'][count-1]-57008.],[m_pk],yerr=[m_err1],fmt='o', color=col)

  ylabel(r'Peak $m_{QSR}$, mag.')

  figure(3)
  subplot(211)
  if(filt=='Sloan_r'):
	errorbar([filename_table['mjd'][count-1]-57008.],[m_pk],yerr=[m_err1],fmt='ro')
  if(filt=='Sloan_g'):
	errorbar([filename_table['mjd'][count-1]-57008.],[m_pk],yerr=[m_err1],fmt='go')
  subplot(212)
  if(filt=='Sloan_r'):
	errorbar([filename_table['mjd'][count-1]-57008.],[pk1],yerr=[sqrt(pk1)],fmt='ro')
  if(filt=='Sloan_g'):
	errorbar([filename_table['mjd'][count-1]-57008.],[pk1],yerr=[sqrt(pk1)],fmt='go')
  if(count != count):
	  figure()
	  imshow(obj.image, cmap='gray', interpolation='None')
	  colorbar()
	  show()
  if(count==60): break

figure(10)
subplot(513)
gca().invert_yaxis()
subplot(514)
gca().invert_yaxis()
subplot(515)
gca().invert_yaxis()
xlabel('Days since 12/17/2014')
subplots_adjust(hspace=0.6, left=0.15, right=0.95, top=0.95, bottom=0.05)


figure(1)
subplot(211)
xlabel('APASS Sloan_r Magnitude')
ylabel(r'$m_{APASS} \ + \  2.5 log_{10}(N_{CCD}/t)$', fontsize=14)
subplot(212)
xlabel('APASS Sloan_g Magnitude')
ylabel(r'$m_{APASS} \ + \  2.5 log_{10}(N_{CCD}/t)$', fontsize=14)


figure(2)
#xlim(0.,2.5)
xlabel('Days since 12/17/2014')
ylabel(r'Zero Point = $m_{APASS} \ + \  2.5 log_{10}(N_{CCD}/t)$', fontsize=14)
grid(True)

figure(3)
subplot(211)
#xlim(0.,2.5)
xlabel('Days since 12/17/2014')
ylabel(r'$Peak m_{QSR}, magnitudes$', fontsize=14)
grid(True)
show()


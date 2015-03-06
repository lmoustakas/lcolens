from AnalysisLibrary import *
from pylab import *
import datetime

rcParams['figure.facecolor']='white'

# SUPER CRUDE PHOTOMETRY

# 1. Read APASS Catalogue
# 2. Loop through data
#    a. Loop through APASS sources
#    b. get the peak value of the image
#    c. compare magnitude to APASS magnitude
#    d. derive correction
#    e. look at mean and rms of correction.
#    f. get quasar image
#    g. get peak values of the four images
#    e. correct the values to magnitudes
#    f. derive time-series with root N errors.

APASS_table = ascii.read('../../data/HE0435_LCOGT/APASS_0438_list.csv')
filename_table = ascii.read('time_ordered_file_list.dat')
dirc = '/data2/romerowo/lcogt_data/he045-1223_wcs_corrected/'
print APASS_table

#ra = 04:38:14.9000 
ra = (4.+38./60.+14.9/60./60.)/24*360 
#dec = -12:17:14.40
dec = -12. - (17./60 +14.4/60./60.)


count = 0
rando = np.random.randint(37,400)
for fnm in filename_table['filename']:
  count += 1
  FM = FITSmanager(dirc+fnm)
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
  if(filename_table['mjd'][count-1]<57008.):
	continue
  
  for k in range(0,len(APASS_table)):
    if(APASS_table[filt][k]!='NA' and k!=21 and k!=9 and k!= 22 ): # entry 21 is on the edge of the field of view.
	obj = SourceImage(FM, APASS_table['radeg'][k], APASS_table['decdeg'][k], 30)
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
  pk1 = np.max(obj.image)
  print '\t', pk1, -2.5*log10(pk1/t),  -2.5*log10(pk1/t) + ZP_mean 
  err1 = sqrt(pk1)
  m_pk = -2.5*log10(pk1/t) + ZP_mean 
  m_err1 = m_pk * 1/err1/log(10.)
  figure(3)
  if(filt=='Sloan_r'):
	errorbar([filename_table['mjd'][count-1]-57008.],[m_pk],yerr=[m_err1],fmt='ro')
  if(filt=='Sloan_g'):
	errorbar([filename_table['mjd'][count-1]-57008.],[m_pk],yerr=[m_err1],fmt='go')
  #if(count == rando):
	  #figure()
	  #imshow(obj.image, cmap='gray', interpolation='None')
	  #show()
  #if(count==10): break
figure(1)
subplot(211)
xlabel('APASS Sloan_r Magnitude')
ylabel(r'$m_{APASS} \ + \  2.5 log_{10}(N_{CCD}/t)$', fontsize=14)
subplot(212)
xlabel('APASS Sloan_g Magnitude')
ylabel(r'$m_{APASS} \ + \  2.5 log_{10}(N_{CCD}/t)$', fontsize=14)

figure(2)
xlim(0.,2.5)
xlabel('Days since 12/17/2014')
ylabel(r'$m_{APASS} \ + \  2.5 log_{10}(N_{CCD}/t)$', fontsize=14)
grid(True)
show()


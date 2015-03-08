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

def APASS_zero_points(FM, APASS_table, APASS_rejects, display=False):
  # KEEP TRACK OF WHICH BAND THE IMAGES ARE AND USE THE CORRESPONDING APASS REFERENCE VALUES.
  filt = 'Sloan_r'
  if( FM.hdulist[0].header['FILTER'] == 'gp'):
    filt = 'Sloan_g'
  print '%d'%(n), fnm
  v1 = []
  v2 = []
  t = FM.get_exposure_time()
  print 't', t

  for k in range(0,len(APASS_table)):
    if(APASS_table[filt][k]!='NA' and k not in APASS_rejects): # entry 21 is on the edge of the field of view.
	obj = SourceImage(FM, APASS_table['radeg'][k], APASS_table['decdeg'][k], 30)
	# GET PEAK
	pk0 = np.max(obj.image)
	# GET BACKGROUND SUBTRACTED PEAK
	pk1 = pk0 - N_bkg
	if(len(obj.image)==30):
  	  #imshow(obj.image)
	  print '%d.%d'%(n,k), fnm, len(obj.image), APASS_table['radeg'][k], APASS_table['decdeg'][k]
	  m_I = -2.5*log10(pk1/t)
	  print '\t',filt, APASS_table[filt][k], m_I
	  v1.append(float(APASS_table[filt][k]))
	  v2.append(m_I)
  v1,v2 = zip(*sorted(zip(v1,v2)))
  ZP = array(v1)-array(v2)
  ZP_mean = sum(ZP)/len(ZP)
  ZP_rms = sqrt(sum(ZP**2)/len(ZP) - ZP_mean**2)
  if(display):
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
  return ZP_mean, ZP_rms

def get_peak_magnitude(obj, N_bkg, sigma_read, ZP_mean, ZP_rms):

  # GET EXPOSURE TIME
  t = obj.FM.get_exposure_time()

  # GET IMAGE COUNTS
  pk0 = np.max(obj.image)

  # GET IMAGE COUNT UNCERTAINTY
  err0 = sqrt(pk0 + (sigma_read)**2)

  # BACKGROUND SUBTRACTED IMAGE COUNTS
  pk1 = pk0 - N_bkg

  # CONVERT TO INSTRUMENTAL MAGNITUDE
  m_pk1 = -2.5*log10(pk1/t)
  err1 = err0/t
  m_err1 =  2.5/err0/log(10.)

  # CONVERT TO CALIBRATED MAGNITUDE
  m_pk2 = -2.5*log10(pk1/t) + ZP_mean 
  m_err2 = sqrt(m_err1**2 + ZP_rms**2 )

  return m_pk2, m_err2

# SET DATA DIRECTORY
dirc = '/data2/romerowo/lcogt_data/he045-1223_wcs_corrected/'

# LOAD APASS SOURCES
APASS_table = ascii.read('../../data/HE0435_LCOGT/APASS_0438_list.csv')

# LOAD TIME ORDERED DATA WITH BACKGROUND COUNTS AND READ NOISE ESTIMATES
filename_table = ascii.read('t_ord_image_stats.dat')
#print APASS_table

# DECIMAL RA and DEC VALUES OF HE0435-1223
#ra = (4.+38./60.+14.9/60./60.)/24*360 
#dec = -12. - (17./60 +14.4/60./60.)

# DECIMAL RA and DEC VALUES OF STAR NEAR HE0435-1223
#ra = (4.+38./60.+14.44/60./60.)/24*360
#dec = -12. - (16./60 +25.4/60./60.)

# DECIMAL RA and DEC VALUES OF STAR NEAR HE0435-1223
#ra = (4.+38./60.+14.60/60./60.)/24*360
#dec = -12. - (16./60 +37.0/60./60.)

# DECIMAL RA and DEC VALUES OF STAR NEAR HE0435-1223
#ra = (4.+38./60.+12.97/60./60.)/24*360
#dec = -12. - (17./60 +51.7/60./60.)

# DECIMAL RA and DEC VALUES OF STAR WITH SIMILAR MAGNITUDE TO HE0435-1223
ra = (4.+38./60.+02.86/60./60.)/24*360
dec = -12. - (16./60 +34.4/60./60.)


# DECIMAL RA and DEC VALUES OF STAR WITH SIMILAR MAGNITUDE TO HE0435-1223
ra = (4.+38./60.+36.49/60./60.)/24*360
dec = -12. - (20./60 +08.2/60./60.)

# KEEP A LIST OF APASS SOURCES THAT DO NOT BEHAVE WELL
APASS_rejects = [9, 21, 22] # 9 and 21 seem to be at the edge of the field of view. 22 seems close to the edge, sometimes we don't catch it.

# LOOP THROUGH FITS FILE IMAGES
for n in range(0,len(filename_table['filename'])):
  
  # GET INFORMATION OF FILENAME AND BACKGROUNDS FROM TABLES
  fnm = filename_table['filename'][n]
  N_bkg = float(filename_table['N_bkg'][n])
  sigma_read = float(filename_table['read_noise'][n])

  # LET'S NOT LOOK AT THE TEST OBSERVATIONS
  mjd_start = 57008.  
  if(filename_table['mjd'][n]< mjd_start):
	continue
  if(n==160): break

  # INITIALIZE THE CUSTOMIZED FITS FILE MANAGER
  FM = FITSmanager(dirc+fnm)

  # GET THE AVERAGE ZERO POINT FROM APASS SOURCES
  figure(1)
  ZP_mean, ZP_rms = APASS_zero_points(FM, APASS_table, APASS_rejects, display=True)

  # GET THE EXPOSURE TIME
  t = FM.get_exposure_time()


  # PLOT MEAN ZERO POINTS OBTAINED FROM APASS SOURCES WITH RMS UNCERTAINTIES
  figure(2)
  print '\t', ZP_mean, ZP_rms
  if(FM.hdulist[0].header['FILTER'] == 'rp'):
	errorbar([filename_table['mjd'][n]-mjd_start],[ZP_mean],yerr=[ZP_rms],fmt='ro')
  if(FM.hdulist[0].header['FILTER'] == 'gp'):
	errorbar([filename_table['mjd'][n]-mjd_start],[ZP_mean],yerr=[ZP_rms],fmt='go')

  # APASS SOURCE CLOSED LOOP TEST

  # PICK A RANDOM APASS SOURCE AND AVOID THE REJECTS
  #random_APASS_index = np.random.randint(0,len(APASS_table['decdeg']))
  try:
	random_APASS_index
  except:
	print 'random APASS index not defined'
  	random_APASS_index = np.random.randint(0,len(APASS_table['decdeg']))
  	while(random_APASS_index in APASS_rejects):
    		random_APASS_index = np.random.randint(0,len(APASS_table['decdeg']))
  obj = SourceImage(FM, APASS_table['radeg'][random_APASS_index], APASS_table['decdeg'][random_APASS_index], 30)

  t_obs = filename_table['mjd'][n]-mjd_start

  figure(3, figsize=(7,12))
  m, m_err = get_peak_magnitude(obj, N_bkg, sigma_read, ZP_mean, ZP_rms)
  # INTERMEDIATE STEPS ADDED HERE FOR DISPLAY
  pk0 = np.max(obj.image)
  err0 = sqrt(pk0 + (sigma_read)**2)
  pk1 = pk0 - N_bkg
  m_pk1 = -2.5*log10(pk1/t)
  err1 = err0/t
  m_err1 =  2.5/err0/log(10.)
  m_pk2 = -2.5*log10(pk1/t) + ZP_mean 
  m_err2 = sqrt(m_err1**2 + ZP_rms**2 )

  col = 'red'
  if(obj.FM.hdulist[0].header['FILTER'] == 'gp'): col='green'
  ax = subplot(611)
  ax.set_yscale('log')
  title('APASS source %d'%random_APASS_index)
  errorbar([t_obs],[N_bkg],yerr=[sqrt(N_bkg+sigma_read**2)],fmt='s', color=col, label='QSR peak')
  errorbar([t_obs],[pk0],yerr=[err0],fmt='*', color=col, label='Background')
  ylabel('Raw Counts')

  subplot(612)
  errorbar([t_obs], [pk1],yerr=[err0],fmt='o', color=col)
  ylabel('$N_{pk}-N_{bkg}$,\nRaw Counts')

  subplot(613)
  errorbar([t_obs], [m_pk1],yerr=[m_err1],fmt='o', color=col)
  ylabel(r'$m_{I}=-2.5 log_{10}((N_{pk}-N_{bkg})/t)$')

  subplot(614)
  errorbar([t_obs], [ZP_mean],yerr=[ZP_rms],fmt='o', color=col)
  ylabel(r'Zero Point')

  if(obj.FM.hdulist[0].header['FILTER'] == 'rp'): 
	  subplot(615)
	  errorbar([t_obs],[m],yerr=[m_err],fmt='o', color=col)
	  val = float(APASS_table['Sloan_r'][random_APASS_index])
	  err_val = float(APASS_table['r_err'][random_APASS_index])
  if(obj.FM.hdulist[0].header['FILTER'] == 'gp'): 
	  subplot(616)
	  errorbar([t_obs],[m],yerr=[m_err],fmt='o', color=col)
	  val = float(APASS_table['Sloan_g'][random_APASS_index])
	  err_val = float(APASS_table['gerr'][random_APASS_index])
  print val, err_val
  errorbar([t_obs], [val], yerr=[err_val], fmt='_', color='gray', linewidth=10, alpha=0.3)
  ylabel(r'Peak mag.')
  gca().invert_yaxis()
  xlabel('Days since mjd %1.0f'%mjd_start)
  subplots_adjust(hspace=0.6, left=0.15, right=0.95, top=0.95, bottom=0.05)


  # GET THE QUASAR IMAGE
  obj = SourceImage(FM, ra, dec, 30)
  m, m_err = get_peak_magnitude(obj, N_bkg, sigma_read, ZP_mean, ZP_rms)

  # INTERMEDIATE STEPS ADDED HERE FOR DISPLAY
  pk0 = np.max(obj.image)
  err0 = sqrt(pk0 + (sigma_read)**2)
  pk1 = pk0 - N_bkg
  m_pk1 = -2.5*log10(pk1/t)
  err1 = err0/t
  m_err1 =  2.5/err0/log(10.)
  m_pk2 = -2.5*log10(pk1/t) + ZP_mean 
  m_err2 = sqrt(m_err1**2 + ZP_rms**2 )

  # PLOT BACKGROUNDS AND QSR IMAGE RAW PEAK COUNTS
  figure(10, figsize=(7,12))
  col = 'red'
  if(obj.FM.hdulist[0].header['FILTER'] == 'gp'): col='green'
  ax1 = subplot(511)
  ax1.set_yscale('log')
  title('HE0434-1223')
  errorbar([t_obs],[N_bkg],yerr=[sqrt(N_bkg+sigma_read**2)],fmt='s', color=col, label='QSR peak')
  errorbar([t_obs],[pk0],yerr=[err0],fmt='*', color=col, label='Background')
  ylabel('Raw Counts')

  subplot(512)
  errorbar([t_obs], [pk1],yerr=[err0],fmt='o', color=col)
  ylabel('$N_{pk}-N_{bkg}$,\nRaw Counts')

  subplot(513)
  errorbar([t_obs], [m_pk1],yerr=[m_err1],fmt='o', color=col)
  ylabel(r'$m_{I}=-2.5 log_{10}((N_{pk}-N_{bkg})/t)$')

  subplot(514)
  errorbar([t_obs], [ZP_mean],yerr=[ZP_rms],fmt='o', color=col)
  ylabel(r'Zero Point')

  subplot(515)
  errorbar([t_obs],[m],yerr=[m_err],fmt='o', color=col)
  ylabel(r'Peak mag.')
  xlabel('Days since mjd %1.0f'%mjd_start)
  subplots_adjust(hspace=0.6, left=0.15, right=0.95, top=0.95, bottom=0.05)
print APASS_table
# FIGURE BEAUTIFICATION


figure(1)
subplot(211)
xlabel('APASS Sloan_r Magnitude')
ylabel(r'$m_{APASS} \ + \  2.5 log_{10}(N_{CCD}/t)$', fontsize=14)
subplot(212)
xlabel('APASS Sloan_g Magnitude')
ylabel(r'$m_{APASS} \ + \  2.5 log_{10}(N_{CCD}/t)$', fontsize=14)

figure(3)
subplot(613)
gca().invert_yaxis()
subplot(614)
gca().invert_yaxis()
subplot(615)
gca().invert_yaxis()
subplot(616)
gca().invert_yaxis()

figure(10)
subplot(513)
gca().invert_yaxis()
subplot(514)
gca().invert_yaxis()
subplot(515)
gca().invert_yaxis()

figure(2)
#xlim(0.,2.5)
xlabel('Days since 12/17/2014')
ylabel(r'Zero Point = $m_{APASS} \ + \  2.5 log_{10}(N_{CCD}/t)$', fontsize=14)
grid(True)

show()


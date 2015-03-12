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

def estimate_total_light(obj, N_bkg, sigma_read, display=False):
  estimate = np.sum(obj.image -  N_bkg)
  uncertainty = np.sqrt(np.sum(obj.image**2)) 
  maskval = np.sqrt((3.*np.sqrt(N_bkg))**2 + sigma_read**2)
  Hmasked = np.ma.masked_where((obj.image-N_bkg)<maskval,obj.image-N_bkg)
  Hmasked2 = np.ma.masked_where((obj.image-N_bkg)<maskval,obj.image)
  estimate = np.sum(Hmasked)
  uncertainty = np.sqrt(np.sum(Hmasked2))
  #if(display or uncertainty/estimate>3e-2):
  if(display):
	  print uncertainty/estimate
	  figure()
	  subplot(221)
          m1 = np.min(obj.image)
          m2 = np.max(obj.image)
          m3 = np.min(obj.image-N_bkg)
          m4 = np.max(obj.image-N_bkg)
	  imshow(obj.image, interpolation='none', vmin=m1, vmax=m2)
	  colorbar()
	  title('Image')
	  subplot(222)
	  imshow(obj.image-N_bkg, interpolation='none', vmin=m3, vmax=m4)
	  colorbar()
	  title('Bkg Subtracted')
	  subplot(223)
	  imshow(Hmasked2, interpolation='none', vmin=m1, vmax=m2)
	  colorbar()
	  subplot(224)
	  imshow(Hmasked, interpolation='none', vmin=m3, vmax=m4)
	  colorbar()
	  show()
  return estimate, uncertainty


def magnitude(value):
  return -2.5*log10(value)

def APASS_zero_points(FM, APASS_table, APASS_rejects, display=False):
  # KEEP TRACK OF WHICH BAND THE IMAGES ARE AND USE THE CORRESPONDING APASS REFERENCE VALUES.
  filt = 'Sloan_r'
  filt_err = 'r_err'
  if( FM.hdulist[0].header['FILTER'] == 'gp'):
    filt = 'Sloan_g'
    filt_err = 'gerr'
  # MAKE LISTST FOR THE CALIBRATED APASS MAGNITUDES AND THE INSTRUMENT MAGNITUDES
  m_APASS = []
  m_I     = []
  m_APASS_unc = []
  m_I_unc     = []
  # LOOP THROUGH APASS OBJECTS TO FILL LISTS. EXCLUDE PRE-ASSESSED APASS REJECTS LIST
  for k in range(0,len(APASS_table)):
    if(APASS_table[filt][k]!='NA' and APASS_table[filt_err][k]!='0' and k not in APASS_rejects): 
	print 'APASS', k
	obj = SourceImage(FM, APASS_table['radeg'][k], APASS_table['decdeg'][k], 40)
	intg, intg_unc = estimate_total_light(obj,N_bkg, sigma_read, display=False)
	m_I.append(magnitude(intg))
	m_APASS.append(float(APASS_table[filt][k]))
	m_I_unc.append(2.5*intg_unc/intg/log10(e))
	m_APASS_unc.append(float(APASS_table[filt_err][k]))
	#print m_I_unc
	#plt.figure()
	#plt.imshow(obj.image, cmap='jet', interpolation='none')
	#plt.colorbar()
	#plt.title('APASS SOURCE %d, magnitude %1.2f\n min N %d, N_bkg %d'%(k,float(APASS_table[filt][k]),np.min(obj.image), N_bkg))
  #plt.show()
  # ESTIMATE ZERO POINTS FOR EACH APASS SOURCE
  ZP = array(m_APASS)-array(m_I)
  print ZP
  # ESTIMATE THE MEAN AND WEIGHTED RMS OF THE ZERO POINTS 
  ZP_mean = sum(ZP)/len(ZP)
  #ZP_rms = sqrt(sum(ZP**2)/len(ZP) - ZP_mean**2)
  weight = 1./np.sqrt(array(m_I_unc)**2+array(m_APASS_unc)**2)
  ZP_rms = sqrt(sum(weight*(ZP-ZP_mean)**2)/sum(weight))
  print ZP_mean, ZP_rms
  #print '**',ZP_mean, ZP_rms
  # WEIGHTED RMS
  if(display):
  	# SORT BY INCREASING MAGNITUDE FOR EASE OF VIEWING
  	m_APASS_ord,m_I_ord = zip(*sorted(zip(m_APASS, m_I)))
  	m_APASS_ord,m_APASS_unc_ord = zip(*sorted(zip(m_APASS, m_APASS_unc)))
  	m_APASS_ord,m_I_unc_ord = zip(*sorted(zip(m_APASS, m_I_unc)))
	if(filt=='Sloan_r'):
		subplot(211)
		errorbar(m_APASS_ord,array(m_APASS_ord)-array(m_I_ord),xerr=m_APASS_unc_ord, yerr=np.sqrt(array(m_I_unc_ord)**2+array(m_APASS_unc_ord)**2) ,fmt='.-')
	if(filt=='Sloan_g'):
		subplot(212)
		errorbar(m_APASS_ord,array(m_APASS_ord)-array(m_I_ord),xerr=m_APASS_unc_ord, yerr=np.sqrt(array(m_I_unc_ord)**2+array(m_APASS_unc_ord)**2),fmt='.-')
  return ZP_mean, ZP_rms


def get_obj_magnitude(obj, N_bkg, sigma_read, ZP_mean, ZP_rms):

  # GET EXPOSURE TIME
  intg, intg_unc = estimate_total_light(obj, N_bkg, sigma_read)

  # CONVERT TO INSTRUMENTAL MAGNITUDE
  m_I     = magnitude(intg)
  m_I_unc = 2.5/intg_unc/log(10.)

  # CONVERT TO CALIBRATED MAGNITUDE
  m     = m_I + ZP_mean 
  m_unc = sqrt(m_I_unc**2 + ZP_rms**2 )

  return m, m_unc


# SET DATA DIRECTORY
dirc = '/data2/romerowo/lcogt_data/he045-1223_wcs_corrected/'

# LOAD APASS SOURCES
APASS_table = ascii.read('../../data/HE0435_LCOGT/APASS_0438_list.csv')

# LOAD TIME ORDERED DATA WITH BACKGROUND COUNTS AND READ NOISE ESTIMATES
filename_table = ascii.read('t_ord_image_stats.dat')
#print APASS_table

# DECIMAL RA and DEC VALUES OF HE0435-1223
ra_qsr = (4.+38./60.+14.9/60./60.)/24*360 
dec_qsr = -12. - (17./60 +14.4/60./60.)

# DECIMAL RA and DEC VALUES OF STAR NEAR HE0435-1223
#ra = (4.+38./60.+14.44/60./60.)/24*360
#dec = -12. - (16./60 +25.4/60./60.)

# DECIMAL RA and DEC VALUES OF STAR NEAR HE0435-1223
#ra = (4.+38./60.+14.60/60./60.)/24*360
#dec = -12. - (16./60 +37.0/60./60.)

# DECIMAL RA and DEC VALUES OF STAR NEAR HE0435-1223
ra = (4.+38./60.+12.97/60./60.)/24*360
dec = -12. - (17./60 +51.7/60./60.)

# DECIMAL RA and DEC VALUES OF STAR WITH SIMILAR MAGNITUDE TO HE0435-1223
ra = (4.+38./60.+02.86/60./60.)/24*360
dec = -12. - (16./60 +34.4/60./60.)

# DECIMAL RA and DEC VALUES OF STAR WITH SIMILAR MAGNITUDE TO HE0435-1223
ra = (4.+38./60.+36.49/60./60.)/24*360
dec = -12. - (20./60 +08.2/60./60.)

# KEEP A LIST OF APASS SOURCES THAT DO NOT BEHAVE WELL
APASS_rejects = [9, 21, 22] # 9 and 21 seem to be at the edge of the field of view. 22 seems close to the edge, sometimes we don't catch it.

def photometry_plot(obj,N_bkg,sigma_read, ZP_mean, ZP_rms):
  m, m_err = get_obj_magnitude(obj, N_bkg, sigma_read, ZP_mean, ZP_rms)
  # INTERMEDIATE STEPS ADDED HERE FOR DISPLAY
  '''
  pk0 = np.max(obj.image)
  pk1 = pk0 - N_bkg
  m_pk1 = -2.5*log10(pk1/t)
  err1 = err0/t
  m_err1 =  2.5/err0/log(10.)
  m_pk2 = -2.5*log10(pk1/t) + ZP_mean 
  m_err2 = sqrt(m_err1**2 + ZP_rms**2 )
  '''
  est,unc = estimate_total_light(obj, N_bkg, sigma_read)
  m_I = -2.5*log10(est)
  m_I_err = 2.5*unc/est/log10(e)

  col = 'red'
  if(obj.FM.hdulist[0].header['FILTER'] == 'gp'): col='green'
  ax = subplot(511)
  ax.set_yscale('log')
  errorbar([t_obs],[N_bkg],yerr=[np.sqrt(N_bkg+sigma_read**2)],fmt='s', color=col, label='Background')
  ylabel('Background, Raw Counts')

  ax1 = subplot(512)
  #ax1.set_yscale('log')
  errorbar([t_obs], [est],yerr=[unc],fmt='o', color=col)
  ylabel('Estimated Illumination from Source,\nRaw Integrated Counts')

  subplot(513)
  errorbar([t_obs], [m_I],yerr=[m_I_err],fmt='o', color=col)
  ylabel(r'$m_{I}$')

  subplot(514)
  errorbar([t_obs], [ZP_mean],yerr=[ZP_rms],fmt='o', color=col)
  ylabel(r'Zero Point Mag.')
  if(obj.FM.hdulist[0].header['FILTER'] == 'rp'): 
	  subplot(515)
	  errorbar([t_obs],[m],yerr=[m_err],fmt='o', color=col)
	  #val = float(APASS_table['Sloan_r'][random_APASS_index])
	  #err_val = float(APASS_table['r_err'][random_APASS_index])
  if(obj.FM.hdulist[0].header['FILTER'] == 'gp'): 
	  subplot(515)
	  errorbar([t_obs],[m],yerr=[m_err],fmt='o', color=col)
	  #val = float(APASS_table['Sloan_g'][random_APASS_index])
	  #err_val = float(APASS_table['gerr'][random_APASS_index])
  #print val, err_val
  #errorbar([t_obs], [val], yerr=[err_val], fmt='_', color='gray', linewidth=10, alpha=0.3)
  ylabel(r'Peak mag.')
  gca().invert_yaxis()
  xlabel('Days since mjd %1.0f'%mjd_start)
  subplots_adjust(hspace=0.6, left=0.15, right=0.95, top=0.95, bottom=0.05)

# LOOP THROUGH FITS FILE IMAGES
count2=0
#for n in range(0,len(filename_table['filename'])):
for n in range(0,len(filename_table['filename'])):
  # GET INFORMATION OF FILENAME AND BACKGROUNDS FROM TABLES
  fnm = filename_table['filename'][n]
  print fnm
  N_bkg = float(filename_table['N_bkg'][n])
  sigma_read = float(filename_table['read_noise'][n])
  print 'Image', n
  # LET'S NOT LOOK AT THE TEST OBSERVATIONS
  mjd_start = 57008.  
  if(filename_table['mjd'][n]< mjd_start):
	continue
  #if(n==60): break
  count2+=1
  # INITIALIZE THE CUSTOMIZED FITS FILE MANAGER
  FM = FITSmanager(dirc+fnm)

  # GET THE AVERAGE ZERO POINT FROM APASS SOURCES
  figure(1)
  ZP_mean, ZP_rms = APASS_zero_points(FM, APASS_table, APASS_rejects, display=True)
  #plt.show()
  # GET THE EXPOSURE TIME
  t = FM.get_exposure_time()

  # PLOT MEAN ZERO POINTS OBTAINED FROM APASS SOURCES WITH RMS UNCERTAINTIES
  figure(2)
  #print '\t', ZP_mean, ZP_rms
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


  t_obs = filename_table['mjd'][n]-mjd_start
  random_APASS_index=5
  obj = SourceImage(FM, APASS_table['radeg'][random_APASS_index], APASS_table['decdeg'][random_APASS_index], 30)  
  figure(3, figsize=(7,12))
  photometry_plot(obj,N_bkg,sigma_read, ZP_mean, ZP_rms)
  # ADD PHOTOMETRIC ESTIMATES FROM APASS SOURCE
  if(obj.FM.hdulist[0].header['FILTER'] == 'rp'): 
	  subplot(515)
	  val = float(APASS_table['Sloan_r'][random_APASS_index])
	  err_val = float(APASS_table['r_err'][random_APASS_index])
  if(obj.FM.hdulist[0].header['FILTER'] == 'gp'): 
	  subplot(515)
	  val = float(APASS_table['Sloan_g'][random_APASS_index])
	  err_val = float(APASS_table['gerr'][random_APASS_index])
  errorbar([t_obs], [val], yerr=[err_val], fmt='_', color='gray', linewidth=10, alpha=0.3)

  ###################################################################################
  obj = SourceImage(FM, ra, dec, 30)  
  figure(4, figsize=(7,12))
  photometry_plot(obj,N_bkg,sigma_read, ZP_mean, ZP_rms)
  '''
  if(count2==9):
	plt.figure()
	plt.imshow(obj.image, interpolation='none')
	plt.colorbar()
	plt.show()
  '''

  # GET THE QUASAR IMAGE
  obj = SourceImage(FM, ra_qsr, dec_qsr, 30)
  intg, intg_unc = estimate_total_light(obj,N_bkg, sigma_read, display=False)

  figure(5, figsize=(7,12))
  photometry_plot(obj,N_bkg,sigma_read, ZP_mean, ZP_rms)
  '''
  if(count2==11):
	plt.figure()
	estimate_total_light(obj, N_bkg, sigma_read, display=True)
	plt.figure()
	levels=arange(0,np.max(obj.image),100)[::-1]
	plt.imshow(obj.image-N_bkg, interpolation='none')
	plt.colorbar()
	plt.show()
  '''
# FIGURE BEAUTIFICATION

figure(1)
subplot(211)
xlabel('APASS Sloan_r Magnitude')
ylabel(r'$m_{APASS} \ + \  2.5 log_{10}(N_{CCD}/t)$', fontsize=14)
subplot(212)
xlabel('APASS Sloan_g Magnitude')
ylabel(r'$m_{APASS} \ + \  2.5 log_{10}(N_{CCD}/t)$', fontsize=14)

figure(3)
subplot(511)
title('APASS source %d'%random_APASS_index)
subplot(513)
gca().invert_yaxis()
subplot(514)
gca().invert_yaxis()
subplot(515)
gca().invert_yaxis()

figure(4)
subplot(511)
title('Star ra: %1.3f dec: %1.3f'%(ra,dec))
subplot(513)
gca().invert_yaxis()
subplot(514)
gca().invert_yaxis()
subplot(515)
gca().invert_yaxis()

figure(5)
subplot(511)
title('HE0435-1223')
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


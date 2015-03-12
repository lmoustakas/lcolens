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

def estimate_total_light(obj, N_bkg):
  estimate = np.sum(obj.image -  N_bkg)
  uncertainty = np.sqrt(np.sum(obj.image**2)) 

  
  #COOKIE CUTTER, WORK IN PROGRESS
  x_max, y_max = np.unravel_index(obj.image.argmax(), obj.image.shape)
  r=[]
  amp=[]
  amp2=[]
  for i in range(0,len(obj.image)):
      for j in range(0,len(obj.image)):
	r.append(np.sqrt((float(i)-x_max)**2 + (float(j)-y_max)**2))
	amp.append(obj.image[i][j]-N_bkg)
	amp2.append(obj.image[i][j]-N_bkg)
  r,amp = zip(*sorted(zip(r, amp)))

  # GET POINT WHERE TOTAL LIGHT CONTRIBUTION IS CHANGING BY < 0.01 %
  # ONE REFINEMENT IS TO MASK VALUES LESS THAN ZERO (OR SQRT(N_BKG)), THEN INTEGRATE RADIALLY OUTWARD
  v = diff(cumsum(amp))/cumsum(amp)[:-1]
  #v_stp_arg = np.argmin((v-1.e-3)**2)
  zero_crossings = np.where(np.diff(np.sign(v-1.e-3)))[0]
  v_stp_arg = zero_crossings[0]
  estimate=cumsum(amp)[v_stp_arg]
  uncertainty=np.sqrt(cumsum(array(amp2)**2)[v_stp_arg])
  #print 'r_stop', r[v_stp_arg]
  
  '''
  #if(np.random.randint(0,100)==0):
  if(2.5*uncertainty/estimate/log10(e) > 0.3):
	  print '*x*', 2.5*uncertainty/estimate/log10(e)
	  plt.figure()
	  plt.imshow(obj.image,interpolation='none',cmap='jet')
	  plt.colorbar()
	  plt.figure()
	  subplot(311)
	  plt.plot(r,amp)
	  subplot(312)
	  plt.semilogy(r[:-1],diff(cumsum(amp))/cumsum(amp)[:-1])
	  print v_stp_arg, v[v_stp_arg]
	  plot([r[v_stp_arg],r[v_stp_arg]],[1.e-7,1.], '--')
	  subplot(313)
	  plt.plot(r,cumsum(amp))
	  plt.show()
  '''

  #print 'e,u', estimate, uncertainty, N_bkg, 2.5*uncertainty/estimate/log10(e)
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
	obj = SourceImage(FM, APASS_table['radeg'][k], APASS_table['decdeg'][k], 40)
	intg, intg_unc = estimate_total_light(obj,N_bkg)
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
  # ESTIMATE THE MEAN AND WEIGHTED RMS OF THE ZERO POINTS 
  ZP_mean = sum(ZP)/len(ZP)
  #ZP_rms = sqrt(sum(ZP**2)/len(ZP) - ZP_mean**2)
  weight = 1./np.sqrt(array(m_I_unc)**2+array(m_APASS_unc)**2)
  ZP_rms = sqrt(sum(weight*(ZP-ZP_mean)**2)/sum(weight))
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
  intg, intg_unc = estimate_total_light(obj, N_bkg)

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
  ax = subplot(511)
  ax.set_yscale('log')
  errorbar([t_obs],[N_bkg],yerr=[sqrt(N_bkg+sigma_read**2)],fmt='s', color=col, label='QSR peak')
  errorbar([t_obs],[pk0],yerr=[err0],fmt='*', color=col, label='Background')
  ylabel('Raw Counts')

  subplot(512)
  errorbar([t_obs], [pk1],yerr=[err0],fmt='o', color=col)
  ylabel('$N_{pk}-N_{bkg}$,\nRaw Counts')

  subplot(513)
  errorbar([t_obs], [m_pk1],yerr=[m_err1],fmt='o', color=col)
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
  N_bkg = float(filename_table['N_bkg'][n])
  sigma_read = float(filename_table['read_noise'][n])
  print 'Image', n
  # LET'S NOT LOOK AT THE TEST OBSERVATIONS
  mjd_start = 57008.  
  if(filename_table['mjd'][n]< mjd_start):
	continue
  if(n==60): break
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
  #if(count2==10):
	#plt.figure()
	#plt.imshow(obj.image)
	#plt.show()
  figure(4, figsize=(7,12))
  photometry_plot(obj,N_bkg,sigma_read, ZP_mean, ZP_rms)

  # GET THE QUASAR IMAGE
  obj = SourceImage(FM, ra_qsr, dec_qsr, 30)
  figure(5, figsize=(7,12))
  photometry_plot(obj,N_bkg,sigma_read, ZP_mean, ZP_rms)
  '''
  if(count2==1):
	plt.figure()
	levels=arange(0,np.max(obj.image-N_bkg),100)[::-1]
	plt.imshow(obj.image-N_bkg, interpolation='none')
	plt.colorbar()
	cs = plt.contour(obj.image-N_bkg, levels=levels, linewidths=0.5, colors='k')
	p = cs.collections[0].get_paths()[0]
	print 'p',p
	v = p.vertices
	print 'v',v
	y = v[:,1]
	#figure()
	#plot(v[:,0],v[:,1])
	print 'y',y
	print 'sum(abs(y))*0.05', sum(abs(y))*0.05
	#s[i,j]=sum(abs(y))*0.05
	#figure()
	#imshow(p)
	#plt.colorbar()
	#print CS
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


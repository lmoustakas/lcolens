from pylab import *
import os
import AnalysisLibrary as AL

rcParams['figure.facecolor']='white'
rcParams['font.size']=14

# read npz file
#figure()
fnm_list = []
# CREATE DATA LISTS FOR APASS AND FOR STARS
# DATA_APPASS = ['APASS_index', 'time_array', 'flux_array', 'flux_uncertainty', 'chisq_array', 'max_chi_array']
APASS_indices = [18,20,23,25,26,27,28,29]
# DATA_APASS = [[18,20,23,25,26,27,28,29],[],[],[],[],[]]
#DATA_APASS = [[18,20,23,25,26,27,28,29], [[[]]*8]*5  ]

#APASS_indices = APASS_indices[::-1]

NSTARS = 14


for k in range(0,500):
  fnm = '/halo_nobackup/lenssim/romerowo/lcolens_outputs/20161001/image_%03d_results.npz'%k
  if(os.path.isfile(fnm)):
    fnm_list.append(fnm)

def get_APASS_time_series(APASS_index, fnm_list):
  time = []
  flux = []
  sig_flux = []
  chisq = []
  maxchi=[]
  for fnm in fnm_list:
    f = np.load(fnm)
    #print f.keys()
    entry = np.where(np.array(APASS_indices) == APASS_index)[0][0]
    cut = np.logical_and(f['APASS_index_list']==APASS_index, f['filter']=='rp')
    if(len(f['APASS_S_CCD'][cut])!=0): 
        time.append(float(f['mjd_obs']))
        flux.append(float(f['APASS_S_CCD'][cut]))
        sig_flux.append(float(f['APASS_sig_S_CCD'][cut]))
        chisq.append(float(f['APASS_chi_sq'][cut]))
        maxchi.append(float(f['APASS_max_chi'][cut]))
    f.close()
  return np.array(time), np.array(flux), np.array(sig_flux), np.array(chisq), np.array(maxchi)

def get_star_time_series(star_index, fnm_list):
  time = []
  flux = []
  sig_flux = []
  chisq = []
  maxchi=[]
  for fnm in fnm_list:
    f = np.load(fnm)
    entry = star_index
    cut = np.logical_and(f['star_index_list']==star_index, f['filter']=='rp')
    if(len(f['star_S_CCD'][cut])!=0):  
        time.append(float(f['mjd_obs']))
        flux.append(float(f['star_S_CCD'][cut]))
        sig_flux.append(float(f['star_S_CCD_unc'][cut]))
        chisq.append(float(f['star_chi_sq'][cut]))
        maxchi.append(float(f['star_max_chi'][cut]))
    f.close()
  return np.array(time), np.array(flux), np.array(sig_flux), np.array(chisq), np.array(maxchi)

print 'PLOTTING APASS FLUXES'
figure(figsize=(8,12))
#ax = subplot(311)
#ax.set_yscale('log')
#color=cm.rainbow(np.linspace(0,1,8))
color=cm.jet(np.linspace(0,1,NSTARS)[::-1])
count=0
APASS_DATA = []
for APASS_index in APASS_indices[::-1]:
  ax1 = subplot(311)
  time, flux, sig_flux, chisq, maxchi = get_APASS_time_series(APASS_index, fnm_list)
  APASS_DATA.append([APASS_index, time, flux, sig_flux, chisq, maxchi])
  m_I = AL.flux2magnitude(flux)
  m_I_unc = AL.fluxErr2magErr(flux, sig_flux)
  #  errorbar(time, flux, yerr=sig_flux, fmt='.')
  #cut = np.logical_and(chisq<1.28, maxchi<4.5)
  cut = np.logical_and(chisq<1.9, maxchi<1.e9)

  errorbar(time[cut]-57008, m_I[cut], yerr=m_I_unc[cut], fmt = ',', color=color[count])
  grid(True)
  ylabel('Instrument Magnitude')
  ax2 = subplot(312, sharex=ax1)
  ax2.set_yscale('log')
  plot(time[cut]-57008, chisq[cut], '.', color=color[count])
  grid(True, which='both')
  ylabel(r'$\chi^2/ndof$')
  ax3 = subplot(313, sharex=ax1)
  ax3.set_yscale('log')
  plot(time[cut]-57008, maxchi[cut], '.', color=color[count])
  grid(True, which='both')
  ylabel(r'max$|\chi|$')
  xlabel('mjd - 57008, days')
  count +=1
ax1 = subplot(311)
ymin, ymax = ax1.get_ylim()
ylim(ymax, ymin)

savefig('APASS_fluxes.png', dpi=50)

print 'PLOTTING STAR FLUXES'
figure(2, figsize=(8,12))
figure(3, figsize=(8,12))
figure(4, figsize=(8,12))
count = 0

ref_time, ref_flux, ref_sig_flux, ref_chisq, ref_maxchi = get_star_time_series(8, fnm_list)
ref_m_I = m_I = AL.flux2magnitude(ref_flux)
ref_cut =  np.logical_and(ref_chisq<1.9, ref_maxchi<1.e9)
ref_m_I = ref_m_I[ref_cut]
ref_time = ref_time[ref_cut]

STAR_DATA = []
for star_index in range(NSTARS)[::-1]:
  figure(2)
  ax1 = subplot(311)
  time, flux, sig_flux, chisq, maxchi = get_star_time_series(star_index, fnm_list)
  STAR_DATA.append([time, flux, sig_flux, chisq, maxchi])
  m_I = AL.flux2magnitude(flux)
  m_I_unc = AL.fluxErr2magErr(flux, sig_flux)
  #  errorbar(time, flux, yerr=sig_flux, fmt='.')
  #cut = np.logical_and(chisq<1.28, maxchi<4.5)
  cut = np.logical_and(chisq<1.9, maxchi<1.e9)
  #print 'check:', star_index, len(time), len(time[cut])

  errorbar(time[cut]-57008, m_I[cut], yerr=m_I_unc[cut], fmt = ',', color = color[count])
  grid(True)
  ylabel('Instrument Magnitude')
  ax2 = subplot(312, sharex=ax1)
  ax2.set_yscale('log')
  plot(time[cut]-57008, chisq[cut], '.', color = color[count])
  grid(True, which='both')
  ylabel(r'$\chi^2/ndof$')
  ax3 = subplot(313, sharex=ax1)
  ax3.set_yscale('log')
  plot(time[cut]-57008, maxchi[cut], '.', color = color[count])
  grid(True, which='both')
  ylabel(r'max$|\chi|$')
  xlabel('mjd - 57008, days')
  figure(3)
  errorbar(time[cut]-57008, m_I[cut] - np.median(m_I[cut]), yerr=m_I_unc[cut], fmt = ',', color = color[count])
  figure(4)
  for i in range(0,len(time)):
    for j in range(0,len(ref_time)):
      if(ref_time[j] == time[i]):
        errorbar(time[i]-57008, m_I[i] - ref_m_I[j], yerr=m_I_unc[i], fmt = ',', color = color[count])
  count +=1
figure(2)
ax1 = subplot(311)
ymin, ymax = ax1.get_ylim()
ylim(ymax, ymin)

print 'len(STAR_DATA)',len(STAR_DATA) 
STAR_DATA=np.array(STAR_DATA)
print STAR_DATA.shape
# CAPTURE ALL THE TIMES IN THE DATA
common_time = []
for i in range(0,len(STAR_DATA)):
    for j in range(0,len(STAR_DATA[i][0])):
        if(STAR_DATA[i][0][j] not in common_time):
            common_time.append(STAR_DATA[i][0][j]) 

common_time = sorted(common_time)
figure()
print len(common_time)
plot(common_time)
savefig('star_fluxes.png', dpi=50)

STAR_DATA_ARRAY = np.zeros( (len(STAR_DATA) , 5, len(common_time) ) )
print STAR_DATA_ARRAY.shape

#show()
# NOW WE FILL IN THE EVENLY SAMPLED NUMPY ARRAY
# WE WILL FILL IN DATA FOR NON-EXISTING TIME POINTS BUT GIVE IT A CHISQ OF 1e16 TO MASK IT.
#print STAR_DATA_ARRAY

for i in range(0,len(STAR_DATA)):
    STAR_DATA_ARRAY[i][0][:] = common_time
    j_shift = 0
    for j in range(0,len(common_time)):
        if(j + j_shift > len(common_time)-1):
            continue
        if(j > len(STAR_DATA[i][0])-1 and j+j_shift < len(common_time)):
            STAR_DATA_ARRAY[i][1][j + j_shift] = 1.
            STAR_DATA_ARRAY[i][2][j + j_shift] = 1.
            STAR_DATA_ARRAY[i][3][j + j_shift] = 1.e20
            STAR_DATA_ARRAY[i][4][j + j_shift] = 1.e20
        if(j < len(STAR_DATA[i][0])):
            while(STAR_DATA_ARRAY[i][0][j + j_shift] != STAR_DATA[i][0][j]):
                STAR_DATA_ARRAY[i][1][j + j_shift] = 1.
                STAR_DATA_ARRAY[i][2][j + j_shift] = 1.
                STAR_DATA_ARRAY[i][3][j + j_shift] = 1.e20
                STAR_DATA_ARRAY[i][4][j + j_shift] = 1.e20
                j_shift +=1
            for k in range(1,5):
                STAR_DATA_ARRAY[i][k][j+j_shift] = STAR_DATA[i][k][j]
        print '-------'
        print i,j, j_shift
        print STAR_DATA_ARRAY[i][0][j+j_shift], STAR_DATA[i][0][j]
        print STAR_DATA_ARRAY[i][1][j+j_shift], STAR_DATA[i][1][j]
        print STAR_DATA_ARRAY[i][2][j+j_shift], STAR_DATA[i][2][j]
        print STAR_DATA_ARRAY[i][3][j+j_shift], STAR_DATA[i][3][j]
        print STAR_DATA_ARRAY[i][4][j+j_shift], STAR_DATA[i][4][j]
        print '-------'
        
print STAR_DATA_ARRAY

figure()
for i in range(0,len(STAR_DATA_ARRAY)):
    cut = np.logical_and(STAR_DATA_ARRAY[i][3]<2., STAR_DATA_ARRAY[i][4]<10.)
    cut = np.logical_and(cut, STAR_DATA_ARRAY[i][2]/STAR_DATA_ARRAY[i][1]<0.3)
    #cut = np.logical_and(cut, STAR_DATA_ARRAY[i][1]-STAR_DATA_ARRAY[i][2]>0.)
    ax = subplot(111)
    ax.set_yscale('log')
    errorbar(STAR_DATA_ARRAY[i][0][cut]-57008, STAR_DATA_ARRAY[i][1][cut], yerr=STAR_DATA_ARRAY[i][2][cut], fmt = ',', color=color[i])

np.savez('STAR_DATA_ARRAY.npz', STAR_DATA_ARRAY=STAR_DATA_ARRAY)

show()

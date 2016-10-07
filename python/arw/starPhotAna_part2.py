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
color=cm.jet(np.linspace(0,1,NSTARS)[::-1])

f = np.load('STAR_DATA_ARRAY.npz')
print f.keys()
figure(figsize=(8,12))
time_matrix = []
mag_matrix = []
err_mag_matrix = []
chisq_matrix = []
maxchi_matrix = []
mask_matrix = []
for k in range(8,len(f['STAR_DATA_ARRAY'])):
  time     = np.array(f['STAR_DATA_ARRAY'][k][0])
  flux     = np.array(f['STAR_DATA_ARRAY'][k][1])
  err_flux = np.array(f['STAR_DATA_ARRAY'][k][2])
  chisq    = np.array(f['STAR_DATA_ARRAY'][k][3])
  maxchi   = np.array(f['STAR_DATA_ARRAY'][k][4])


  mag = AL.flux2magnitude(flux)
  #err_mag = AL.fluxErr2magErr(flux, err_flux)
  err_mag = np.abs(err_flux) / flux

  cut1 = chisq<3.
  cut2 = maxchi<5.
  cut3 = flux>0.
  cut4 = flux-err_flux > 0.
  cut5 = np.isfinite(mag)
  cut6 = np.isfinite(err_mag)
  cut = np.logical_and(cut1, cut2)
  cut = np.logical_and(cut, cut3) 
  cut = np.logical_and(cut, cut4)
  cut = np.logical_and(cut, cut5)
  cut = np.logical_and(cut, cut6)


  time_matrix.append(time)
  mag_matrix.append(mag)
  err_mag_matrix.append(err_mag)
  chisq_matrix.append(chisq)
  maxchi_matrix.append(maxchi)
  mask_matrix.append(cut)


  #cut = np.logical_and(chisq<1.e20, isfinite(mag))

  N0 = len(time)
  time = time[cut]-57008.
  flux = flux[cut]
  err_flux = err_flux[cut]
  chisq = chisq[cut]
  maxchi = maxchi[cut]
  mag = mag[cut]
  err_mag = err_mag[cut]

  print 'Star',k,'N0;', N0, 'N', len(time)
  ax1 = subplot(411)
  errorbar(time, mag, yerr=err_mag, fmt='.', lw=1, ms=10, color=color[k])
  ylabel('Instrument Mag')
 
  ax11 = subplot(412)
  errorbar(time, mag - np.mean(mag), yerr=err_mag, fmt='.', lw=1, ms=5, color=color[k])
  ylabel('Instrument Mag - Mean')
  #plot(time, mag-np.mean(mag), '.', lw=2, ms=5, color=color[k])
 
  ax2 = subplot(413, sharex=ax1)
  ax2.set_yscale('log')
  plot(time, chisq, '.-', lw=1, ms=10, color=color[k])
  ylabel('$\chi^2/ndof$')
  ax3 = subplot(414, sharex=ax1)
  ax3.set_yscale('log')
  plot(time, maxchi, '.-', lw=1, ms=10,color = color[k])
  ylabel('max$\chi$')
  xlabel('time, mjd-57008')

#show()
time_matrix    = np.array(time_matrix)
mag_matrix     = np.array(mag_matrix)
err_mag_matrix = np.array(err_mag_matrix)
chisq_matrix   = np.array(chisq_matrix)
maxchi_matrix  = np.array(maxchi_matrix)
mask_matrix = np.array(mask_matrix)


print '**'
print time_matrix.shape
print mag_matrix.shape
print err_mag_matrix.shape

mask_matrix = np.logical_and(np.isfinite(mag_matrix), np.isfinite(err_mag_matrix))
mask_matrix = np.logical_and(mask_matrix, chisq_matrix<3.)
mask_matrix = np.logical_and(mask_matrix, maxchi_matrix<5.)
print 'mask_matrix.shape', mask_matrix.shape 

# REMOVE SLICES WITH NO SUCCESSFUL ESTIMATES OF MAGNITUDE OR MAGNITUDE UNCERTAINTY
print 'np.sum(mask_matrix)', np.sum(mask_matrix)
print 'np.sum(mask_matrix, axis=0)', np.sum(mask_matrix, axis=0)
print 'np.sum(mask_matrix, axis=1)', np.sum(mask_matrix, axis=1)

sum_valid_array = np.sum(mask_matrix, axis=0)

clean_flag=True
while clean_flag:
  print mask_matrix.shape
  sum_valid_array = np.sum(mask_matrix, axis=0)
  for k in range(0, len(sum_valid_array)):
    clean_flag=False
    #print k, sum_valid_array[k]
    if sum_valid_array[k] < 3:
      print k, sum_valid_array[k]
      time_matrix    = np.delete(time_matrix, k, axis=1)
      mag_matrix     = np.delete(mag_matrix, k, axis=1)
      err_mag_matrix = np.delete(err_mag_matrix, k, axis=1)
      maxchi_matrix  = np.delete(maxchi_matrix, k, axis=1)
      mask_matrix    = np.delete(mask_matrix, k, axis=1)
      clean_flag=True
      break
  print 'continuing here', clean_flag


print mask_matrix.shape

figure()
subplot(211)
plot(time_matrix, mag_matrix, '.')
subplot(212)
plot(time_matrix[mask_matrix], mag_matrix[mask_matrix], '.')


I2d, T2d = np.indices(mag_matrix.shape)
print T2d
print I2d
print T2d.shape
print I2d.shape

def func((i2d,t2d), *parms):
  N = len(t2d)
  M = len(i2d)
  #print 'N,M', N,M
  #print t2d
  #print i2d
  #print 'N,M', N,M
  star_mags = np.array(parms[0:np.max(i2d)+1])
  time_var = np.array(parms[np.max(i2d)+1:])
  _mags = np.zeros(N)
  #print 'len(_mags), len(star_mags), len(time_var)', len(_mags), len(star_mags), len(time_var)
  for k in range(0,N):
   #print k, i2d[k], t2d[k]
   _mags[k] = star_mags[i2d[k]] + time_var[t2d[k]]
  return _mags.ravel()

star_mags = np.mean(np.ma.masked_where(mask_matrix, mag_matrix), axis=1)
time_var = np.mean( (np.ma.masked_where(mask_matrix, mag_matrix).transpose() -star_mags).transpose(), axis=0 )

N,M = time_matrix.shape
p0 = np.zeros(N+M)

print 'p0', p0

from scipy.optimize import curve_fit

print 'Check_finite'
print '\ttime_matrix', np.any(np.isneginf(time_matrix))
print '\tmag_matrix.ravel()', np.any(np.isneginf(mag_matrix[mask_matrix]))
print '\terr_mag_matrix.ravel()', np.any(np.isneginf(err_mag_matrix[mask_matrix]))
print 'START CURVE_FIT'
#popt, pcov = curve_fit(func, (T2d,I2d), mag_matrix.ravel(), sigma=err_mag_matrix.ravel(), p0=p0, check_finite=False, absolute_sigma=True)
popt, pcov = curve_fit(func, (I2d[mask_matrix], T2d[mask_matrix]), mag_matrix[mask_matrix], sigma=err_mag_matrix[mask_matrix], p0=p0, check_finite=True, absolute_sigma=True)
#show()
print popt

figure()
plot(popt[N:], '.')

figure()
for k in range(0,len(time_matrix)):
    time    = time_matrix[k] - 57008.
    mag     = mag_matrix[k]
    err_mag = err_mag_matrix[k]  
    mask    = mask_matrix[k]
    subplot(len(time_matrix)+1, 1, k+1)
    errorbar(time[mask], mag[mask] - popt[N:][mask] - popt[k], yerr=err_mag[mask], fmt='.')
    ylim(-0.15, 0.15)
    grid(True)

figure()
for k in range(0,len(time_matrix)):
    time    = time_matrix[k] - 57008.
    mag     = mag_matrix[k]
    err_mag = err_mag_matrix[k]
    mask    = mask_matrix[k]
    subplot(len(time_matrix)+1, 1, k+1)
    #hist(mag[mask])
    hist(err_mag[mask], range=(0,0.1), bins=50)
    xlim(0., 0.1)
    grid(True)

figure()
for k in range(0,len(time_matrix)):
    time    = time_matrix[k] - 57008.
    mag     = mag_matrix[k]
    err_mag = err_mag_matrix[k]
    mask    = mask_matrix[k]
    subplot(len(time_matrix)+1, 1, k+1)
    #hist(mag[mask])
    hist(mag[mask]  - popt[N:][mask] - popt[k], range=(-0.15,0.15), bins=51)
    xlim(-0.15, 0.15)
    grid(True)
show()
show()

print 'END CURVE_FIT'
#print 'popt', popt
#model_matrix = func((T2d,I2d), *popt).reshape(time_matrix.shape)
figure()
for k in range(0,len(time_matrix)):
  subplot(len(time_matrix), 1, k+1)
  errorbar(time_matrix[k]-57008., model_matrix[k] - mag_matrix[k], yerr=err_mag_matrix[k], fmt='.', color=color[k+8])
  ylim(-0.15, 0.15)
  grid(True)
suptitle('fit residuals')
print 'star_mags', star_mags
#figure()
#print len(popt[NSTARS:]), len(np.sqrt(np.diag(pcov)[NSTARS:]))
#errorbar(time_matrix[k]-57008., popt[NSTARS:], yerr=np.sqrt(np.diag(pcov)[NSTARS:]))
figure()
plot(time_matrix, mag_matrix,'.')
figure()
vals = (mag_matrix.transpose() -star_mags).transpose()
print 'vals.shape', vals.shape
#print '*', np.mean(vals, axis=0).shape, np.mean(vals, axis=0)
plot(time_matrix, vals, '.')
#print mag_matrix.ravel()
print 'time_matrix.shape', time_matrix.shape
plot(time_matrix[0],time_var, 'o-', mfc='none', mec='purple', mew=2, ms=10)

show()

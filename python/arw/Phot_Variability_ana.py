#!/usr/bin/env python

'''
2016 October 28
Andrew Romero-Wolf
Post-process photometric estimation
'''

import os 
import AnalysisLibrary as AL
from pylab import *
from astropy.io import ascii

rcParams['figure.facecolor'] = 'white'
rcParams['font.size'] = 16


# Set cuts
num_good_PSF_fit_cut = 6

#Star_chi_max_cut = 4.5
#Star_chi_sq_cut = 1.15
Star_chi_max_cut = 4.3
Star_chi_sq_cut = 1.3
num_Star_good_cut = 30
Star_fractional_uncertainty_cut = 0.1
#Star_fractional_uncertainty_cut = 0.1

# Time variabiluity Cuts
# NUmber of images a star participates in: 600
# Outliers MAD: 2. 

QSR_chi_max_cut = 4.5
QSR_chi_sq_cut = 1.3
num_QSR_good_cut = 1
QSR_fractional_uncertainty_cut = 0.3



###########################################

def filter_cuts(fnm_list, filter = 'rp'):
  new_fnm_list = []
  for fnm in fnm_list:
    f = np.load(fnm)
    if f['filter'] == filter:
        new_fnm_list.append(fnm)
    f.close()
  return(new_fnm_list)

###########################################

def PSF_quality_cuts(fnm_list):
  num_good_PSF = []
  new_fnm_list = []
  figure(figsize=(10,12))
  ax1=subplot(211)
  for fnm in fnm_list:
    f = np.load(fnm)
    num_good_PSF.append(f['num_good_PSF_fit'])
    num = int(fnm.split('image_')[1].split('_')[0])
    plot([num], num_good_PSF[-1], 'ko')
    if(f['num_good_PSF_fit'] >= num_good_PSF_fit_cut):
        new_fnm_list.append(fnm)
    f.close()
  xlabel('Image Number')
  ylabel('Number of PSF Star Fits Passing Quality Cuts')
  x1,x2 = ax1.get_xlim()
  plot([x1,x2],[num_good_PSF_fit_cut, num_good_PSF_fit_cut], 'r--', lw=2)

  grid(True)
  ax2=subplot(212)
  #plot(sorted(num_good_PSF), np.linspace(0,1., len(num_good_PSF)),'ko')
  hist(num_good_PSF, bins=np.linspace(0,13, 14), facecolor='none')
  xlabel('Number of PSF Star Fits Passing Quality Cuts')
  ylabel('Number of Images')
  xlim(-0.1, 13.)
  grid(True)
  y1,y2 = ax2.get_ylim()
  plot([num_good_PSF_fit_cut, num_good_PSF_fit_cut],[y1,y2], 'r--', lw=2)
  return(new_fnm_list)
    
#############################################################################################

def Star_quality_cuts(fnm_list):
  chis = []
  max_chis = []
  vals = []
  new_fnm_list = []
  frac_uncertainty = []
  figure(figsize=(10,12))
  ax1 = subplot(211)
  for fnm in fnm_list:
    f = np.load(fnm)
    chi_list      = f['chi_sq']
    max_chi_list  = f['max_chi']
    S_CCD         = f['S_CCD']
    sig_S_CCD     = f['sig_S_CCD']
    
    [frac_uncertainty.append(fr) for fr in sig_S_CCD/S_CCD]
    [chis.append(c) for c in chi_list ]
    [max_chis.append(m) for m in max_chi_list ]

    qual_cut = np.logical_and(chi_list<Star_chi_sq_cut, max_chi_list<Star_chi_max_cut)
    qual_cut = np.logical_and(sig_S_CCD/S_CCD<Star_fractional_uncertainty_cut, qual_cut)
    qual_cut = np.logical_and(S_CCD>0., qual_cut)
    qual_cut = np.logical_and(sig_S_CCD>0., qual_cut)
    qual_cut = np.logical_and(np.isfinite(S_CCD), qual_cut)
    qual_cut = np.logical_and(np.isfinite(sig_S_CCD), qual_cut)
    vals.append(np.sum(qual_cut))
    if(np.sum(qual_cut)>=num_Star_good_cut):
      new_fnm_list.append(fnm)
    f.close()
    num = int(fnm.split('image_')[1].split('_')[0])
    plot([num], [np.sum(qual_cut)], 'ko')
  xlabel('Image Number')
  ylabel('Number of Star Fits Passing Quality Cuts')
  x1,x2 = ax1.get_xlim()
  plot([x1,x2],[num_Star_good_cut, num_Star_good_cut], 'r--', lw=2)
  grid(True)
  ax2=subplot(212)
  #plot(sorted(num_good_PSF), np.linspace(0,1., len(num_good_PSF)),'ko')
  hist(vals, bins=np.linspace(0,100., 21), facecolor='none')
  xlabel('Number of Star Fits Passing Quality Cuts')
  ylabel('Number of Images')
  xlim(-0.1, np.max(vals)+1)
  grid(True)
  y1,y2 = ax2.get_ylim()
  plot([num_Star_good_cut, num_Star_good_cut],[y1,y2], 'r--', lw=2)
  
  chis = np.array(chis)
  max_chis = np.array(max_chis)

  display_cut = np.isfinite(chis)
  display_cut = np.logical_and( np.isfinite(frac_uncertainty), display_cut )
  display_cut = np.logical_and( frac_uncertainty>0., display_cut )
  #chis = np.array(chis)
  #print 'chis.shape', chis.shape
  #print 'qual_cut.shape',qual_cut.shape
  #chis = chis[qual_cut]

  figure()
  frac_uncertainty = np.array(frac_uncertainty)
  ax = subplot(111)
  ax.set_yscale('log')
  ax.set_xscale('log')
  h,b = np.histogram( np.log10(frac_uncertainty[display_cut]), bins=np.linspace(-3,1.,200))
  xlabel('Star Fractional Uncertainty on Flux')
  title('Star Flux Uncertainties')
  plot(10**b[1:], h, drawstyle='steps')
  y11 = 0.5
  y12 = 2.*np.max(h)
  ylim(y11, y12)
  plot([Star_fractional_uncertainty_cut, Star_fractional_uncertainty_cut],[0.5, 2.*np.max(h)], 'r--', lw=2)
  grid(True, which='both')
  

  figure(figsize=(10,8))
  suptitle('Star Fit $\chi^2$ and max$|\chi|$')
  h,b = np.histogram( np.log10(chis[display_cut]), bins=500)
  ax1 = subplot(221)
  ax1.set_yscale('log')
  ax1.set_xscale('log')
  plot(10**b[1:], h+1.e-3, drawstyle='steps')
  y11 = 0.5
  y12 = 2.*np.max(h)
  ylim(y11, y12)
  x11, x12 = ax1.get_xlim()
  x11 = 10**np.min(b)
  xlim(x11, x12)
  plot([Star_chi_sq_cut, Star_chi_sq_cut],[0.5, 2.*np.max(h)], 'r--', lw=2)
  grid(True, which='both')

  h,b = np.histogram( np.log10(max_chis[display_cut]), bins=500)
  ax2 = subplot(224)
  ax2.set_yscale('log')
  ax2.set_xscale('log')
  plot(10**b[1:], h+1.e-3, drawstyle='steps')
  ylim(0.5, 2.*np.max(h))
  plot([Star_chi_max_cut, Star_chi_max_cut],[0.5, 2.*np.max(h)], 'r--', lw=2)
  x21, x22 = ax2.get_xlim()
  x21 = 2.
  xlim(x21, x22)
  grid(True, which='both')
  xlabel(r'max$|\chi|$')

  ax3 = subplot(223, sharex=ax1)
  nbins = 100
  H, xedges, yedges = np.histogram2d(np.log10(chis[display_cut]), np.log10(max_chis[display_cut]), bins=nbins)
  H = np.rot90(H)
  H = np.flipud(H) 
  Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
 
  from matplotlib.colors import LogNorm
  ax3.pcolormesh(10**xedges,10**yedges,Hmasked, cmap='viridis', norm=LogNorm()) 
  xlabel(r'$\chi^2$')
  ylabel(r'max$|\chi|$')
  xlim(x11,x12)
  ylim(x21,x22)
  plot([Star_chi_sq_cut, Star_chi_sq_cut], [x21, x22], 'r--', lw=2)
  plot([x11, x12], [Star_chi_max_cut, Star_chi_max_cut], 'r--', lw=2)

  ax3.set_yscale('log')
  ax3.set_xscale('log')

  figure()
  #print 'np.min(vals)', np.min(vals)
  print 'len(vals)',len(vals)
  h, b = np.histogram(vals, bins=np.linspace(0,np.max(vals)*1.1, int(1.1*np.max(vals))+1))
  plot(b[1:], np.cumsum(h), drawstyle='steps')
  #plot(b[1:], h, drawstyle='steps')
  grid(True)
  #xlim(-0.5, 6.5)
  ylim(0., 1.1*np.sum(h))
  xlim(0,np.max(vals)*1.1)
  #ylim(0., 1.1*np.max(h))
  plot([num_Star_good_cut-0.5, num_Star_good_cut-0.5],[0., 1.1*np.sum(h)], 'r--', lw=2)
  xlabel('Number of Stars Passing Quality Cut')
  ylabel('Cumulative Counts')
  title('Star Quality Cuts')

  return(new_fnm_list)

######################################################################################################


def QSR_quality_cuts(fnm_list, kind='wg'):
  # kind is either 'wg' or 'ng'
  chis = []
  max_chis = []
  vals = []
  new_fnm_list = []
  frac_uncertainty = []
  for fnm in fnm_list:
    f = np.load(fnm)
    #print f.keys()
    #print f['qsr_wg_parms'][10:14]
    #print np.sqrt(np.diagonal(f['qsr_wg_covar']))[10:14]
    #print np.sqrt(np.diagonal(f['qsr_wg_covar']))[10:14] / f['qsr_wg_parms'][10:14]
    chis.append(f['qsr_%s_chisq'%kind])
    max_chis.append(f['qsr_%s_max_chi'%kind])
    [frac_uncertainty.append(fr) for fr in np.sqrt(np.diagonal(f['qsr_%s_covar'%kind]))[10:14] / f['qsr_%s_parms'%kind][10:14]]

    qual_cut = np.logical_and(f['qsr_%s_chisq'%kind]<QSR_chi_sq_cut, f['qsr_%s_max_chi'%kind]<QSR_chi_max_cut)
    qual_cut = np.logical_and(np.all(np.sqrt(np.diagonal(f['qsr_%s_covar'%kind]))[10:14] / f['qsr_%s_parms'%kind][10:14]<QSR_fractional_uncertainty_cut), qual_cut)
    qual_cut = np.logical_and(np.all(f['qsr_%s_parms'%kind][10:14]>0.), qual_cut)
    qual_cut = np.logical_and(np.all(np.sqrt(np.diagonal(f['qsr_%s_covar'%kind]))[10:14]>0.), qual_cut)
    qual_cut = np.logical_and(np.all(np.isfinite(f['qsr_%s_parms'%kind][10:14])), qual_cut)
    qual_cut = np.logical_and(np.all(np.isfinite(np.sqrt(np.diagonal(f['qsr_%s_covar'%kind]))[10:14])), qual_cut)
    vals.append(np.sum(qual_cut))
    if(np.sum(qual_cut)>=1):
      new_fnm_list.append(fnm)
    f.close()

  chis = np.array(chis)
  max_chis = np.array(max_chis)
  frac_uncertainty = np.array(frac_uncertainty)

  chi_cut = np.isfinite(chis)
  frac_unc_cut = np.logical_and( np.isfinite(frac_uncertainty), frac_uncertainty>0. )

  figure()
  ax = subplot(111)
  ax.set_yscale('log')
  ax.set_xscale('log')
  h,b = np.histogram( np.log10(frac_uncertainty[frac_unc_cut]), bins=np.linspace(-3,1.,100))
  xlabel('QSR Fractional Uncertainty on Flux')
  title('QSR Flux Uncertainties')

  plot(10**b[1:], h, drawstyle='steps')
  y11 = 0.5
  y12 = 2.*np.max(h)
  ylim(y11, y12)
  plot([QSR_fractional_uncertainty_cut, QSR_fractional_uncertainty_cut],[0.5, 2.*np.max(h)], 'r--', lw=2)
  grid(True, which='both')

  figure()
  suptitle('QSR Fit $\chi^2$ and max$|\chi|$')
  h,b = np.histogram( np.log10(chis[chi_cut]), bins=200)
  ax1 = subplot(221)
  ax1.set_yscale('log')
  ax1.set_xscale('log')
  plot(10**b[1:], h+1.e-3, drawstyle='steps')
  y11 = 0.5
  y12 = 2.*np.max(h)
  ylim(y11, y12)
  x11, x12 = ax1.get_xlim()
  x11 = 10**np.min(b)
  xlim(x11, x12)
  plot([QSR_chi_sq_cut, QSR_chi_sq_cut],[0.5, 2.*np.max(h)], 'r--', lw=2)
  grid(True, which='both')

  h,b = np.histogram( np.log10(max_chis[chi_cut]), bins=200)
  ax2 = subplot(224)
  ax2.set_yscale('log')
  ax2.set_xscale('log')
  plot(10**b[1:], h+1.e-3, drawstyle='steps')
  ylim(0.5, 2.*np.max(h))
  plot([QSR_chi_max_cut, QSR_chi_max_cut],[0.5, 2.*np.max(h)], 'r--', lw=2)
  x21, x22 = ax2.get_xlim()
  x21 = 2.
  xlim(x21, x22)
  grid(True, which='both')
  xlabel(r'max$|\chi|$')

  ax3 = subplot(223, sharex=ax1)
  nbins = 100
  H, xedges, yedges = np.histogram2d(np.log10(chis[chi_cut]), np.log10(max_chis[chi_cut]), bins=nbins)
  H = np.rot90(H)
  H = np.flipud(H) 
  Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
 
  from matplotlib.colors import LogNorm
  ax3.pcolormesh(10**xedges,10**yedges,Hmasked, cmap='viridis', norm=LogNorm()) 
  xlabel(r'$\chi^2$')
  ylabel(r'max$|\chi|$')
  xlim(x11,x12)
  ylim(x21,x22)
  plot([QSR_chi_sq_cut, QSR_chi_sq_cut], [0.1, 1000.], 'r--', lw=2)
  plot([0.1, 1000.], [QSR_chi_max_cut, QSR_chi_max_cut], 'r--', lw=2)

  ax3.set_yscale('log')
  ax3.set_xscale('log')

  figure()
  #print 'np.min(vals)', np.min(vals)
  h, b = np.histogram(vals, bins=np.linspace(-1.5,1.5, 4))
  plot(b[1:], np.cumsum(h), drawstyle='steps')
  grid(True)
  xlim(-0.5, 1.5)
  ylim(0., 1.1*np.sum(h))
  plot([num_QSR_good_cut-0.5, num_QSR_good_cut-0.5],[0., 1.1*np.sum(h)], 'r--', lw=2)
  xlabel('Number of QSR Passing Quality Cut')
  ylabel('Cumulative Counts')
  title('QSR Quality Cuts')


  return(new_fnm_list)

######################################################################################################

def check_errors(err_fnm):
    lc = 0
    for line in file(err_fnm):
        lc+=1
    return lc

######################################################################################################

def binned_light_curve(_times, _mags, _mag_err, bins):
  # print 'bins', bins
  t_binned = []
  wmean_binned = []
  werr_binned = []
  for i in range(0,len(bins)-1):
    cut = np.logical_and(_times>=bins[i], _times<bins[i+1])
    N = np.sum(cut)
    #print 'N', N, i, bins[i], bins[i+1], np.min(_times), np.max(_times)
    if(N!=0):
        #print '\tN =',N
        t_binned.append((bins[i+1] + bins[i])/2.)
        #print i, N, bins[i]
        #print '\t', np.min(times[cut]), np.max(times[cut])
        #ax1 = subplot(411)
        weights = 1./(_mag_err[cut]**2)
        wmean = np.sum( weights * _mags[cut] ) / np.sum( weights )
        werr = np.sqrt( 1. / np.sum( weights ) )
        #print wmean, werr
        wmean_binned.append(wmean)
        werr_binned.append(werr)
        #print '\tlen(wmean_binned) =',len(wmean_binned)
        #errorbar((bins[i+1] + bins[i])/2., wmean_1, yerr=werr_1, fmt='b.')
  # print len(wmean_binned)
  return np.array(t_binned), np.array(wmean_binned), np.array(werr_binned)

######################################################################################################


def flux_corrected_curves(fnm_list, t0=57008.):
  star_table = ascii.read('sample.csv')
  #print star_table
  #print star_table['SDSS_r']

  print 'Function: flux_corrected_curves'
  # index the times, magnitudes, errors, and create mask matrix for the values that pass our cuts.
  times = [] 
  image_num = []
  f = np.load(fnm_list[0])
  N_stars = len(f['S_CCD'])
  f.close()
  star_mag_matrix = np.zeros((N_stars,len(fnm_list)))
  star_mag_err_matrix = np.zeros((N_stars,len(fnm_list)))
  star_mask_matrix = np.zeros((N_stars,len(fnm_list)), dtype = bool)
  #for fnm in fnm_list:
  for i in range(0,len(fnm_list)):
    f = np.load(fnm_list[i])
    #print fnm_list[i]
    times.append(f['mjd_obs']-t0)
    num = int(fnm_list[i].split('image_')[1].split('_')[0])
    image_num.append(num)
    S_CCD     = f['S_CCD']
    sig_S_CCD = f['sig_S_CCD']
    chi_sq    = f['chi_sq']
    max_chi   = f['max_chi']
#AL.flux2magnitude(qsr_mag_matrix[k]) - popt[N:], yerr=AL.fluxErr2magErr(qsr_mag_matrix[k],
    for j in range(0,len(S_CCD)):
      star_mag_matrix[j][i]     = AL.flux2magnitude( S_CCD[j] )
      star_mag_err_matrix[j][i] = AL.fluxErr2magErr( S_CCD[j], sig_S_CCD[j] )
      fr = sig_S_CCD[j]/ S_CCD[j]
      if(      chi_sq[j]<Star_chi_sq_cut
           and max_chi[j]<Star_chi_max_cut
           and fr < Star_fractional_uncertainty_cut
           and S_CCD[j]>0.):
             star_mask_matrix[j][i] = True
    f.close()


  ########################################################################################
  # PURGE STARS THAT WON'T CONTRIBUTE TO THE FIT
  num_image_contributions_by_star = np.sum(star_mask_matrix, axis=1)
  star_cut = num_image_contributions_by_star > 600
  star_index = np.where(star_cut)[0]
  
  figure(figsize=(10,12))
  ax1 = subplot(211)
  plot(num_image_contributions_by_star, 'ko')
  plot(np.arange(0,len(num_image_contributions_by_star)), 600*np.ones(len(num_image_contributions_by_star)), 'r--', lw=2)
  xlabel('Star Index')
  ylabel('Number of Images\nPassing Quality Cuts')
  grid(True)
  title('Star Time Series Cuts')
  xlim(0, len(num_image_contributions_by_star))
  ax2 = subplot(212)
  hist(num_image_contributions_by_star, bins=np.linspace(0.,900,46), cumulative=-1, facecolor='none')
  y1, y2 = ax2.get_ylim()
  plot([600, 600], [y1,y2], 'r--', lw=2)
  xlabel('Number of Images Passing Quality Cuts')
  ylabel('Cumulative Counts')
  grid(True)

  ########################################################################################

  N_stars = np.sum(star_cut)
  print '\n\tNum. Images Where a star contributes > 600 Times:', N_stars
  print '\n\tReshape Star Data Matrix from', star_mag_matrix.shape, 'to', star_mag_matrix[star_cut,:].shape

  star_mag_matrix     = star_mag_matrix[star_cut,:]
  star_mag_err_matrix = star_mag_err_matrix[star_cut,:]
  star_mask_matrix    = star_mask_matrix[star_cut,:]

  times = np.array(times)
  image_num = np.array(image_num)

  ########################################################################################

  print '\n\tMaking Initial Estimates of Average Fluxes and Time Variability'
  star_mags = np.zeros(N_stars)
  for ii in range(0,N_stars):
    star_mags[ii] = np.median(star_mag_matrix[ii,:][star_mask_matrix[ii,:]])
    #star_mags[ii] = np.mean(star_mag_matrix[ii,:][star_mask_matrix[ii,:]])

  time_var = np.zeros(star_mag_matrix.shape[1])
  for jj in range(0, star_mag_matrix.shape[1]):
    res = star_mag_matrix[:,jj] - star_mags[:]
    time_var[jj] = np.median(res[star_mask_matrix[:,jj]])
    #time_var[jj] = np.mean(res[star_mask_matrix[:,jj]])

  print '\n\tCheck Magnitudes and Time Series were all Estimates'
  print '\t\tnp.sum(star_mags==0.)',np.sum(star_mags==0.)
  print '\t\tnp.sum(time_var==0.)',np.sum(time_var==0.)

 
  ########################################################################################
  # REMOVE TIME POINTS WHERE NO STARS PARTICIPATE

  time_cut = time_var!=0.
  print '\n\tReshape Star Data Matrix from', star_mag_matrix.shape, 'to', star_mag_matrix[:,time_cut].shape

  time_index = np.where(time_cut)[0]
  star_mag_matrix     = star_mag_matrix[:,time_cut]
  star_mag_err_matrix = star_mag_err_matrix[:,time_cut]
  star_mask_matrix    = star_mask_matrix[:,time_cut]
  time_var            = time_var[time_cut]
  times               = times[time_cut]
  image_num           = image_num[time_cut]

  ########################################################################################

  print '\n\tRe-making Initial Estimates of Average Fluxes and Time Variability'
  star_mags = np.zeros(N_stars)
  for ii in range(0,N_stars):
    star_mags[ii] = np.median(star_mag_matrix[ii,:][star_mask_matrix[ii,:]])
    #star_mags[ii] = np.mean(star_mag_matrix[ii,:][star_mask_matrix[ii,:]])

  time_var = np.zeros(star_mag_matrix.shape[1])
  for jj in range(0, star_mag_matrix.shape[1]):
    res = star_mag_matrix[:,jj] - star_mags[:]
    time_var[jj] = np.median(res[star_mask_matrix[:,jj]])
    #time_var[jj] = np.mean(res[star_mask_matrix[:,jj]])

  print '\n\tCheck Magnitudes and Time Series were all Estimates'
  print '\t\tnp.sum(star_mags==0.)',np.sum(star_mags==0.)
  print '\t\tnp.sum(time_var==0.)',np.sum(time_var==0.)

  '''
  num_image_contributions_by_star = np.sum(star_mask_matrix, axis=1)
  star_cut = num_image_contributions_by_star > 600
  star_index = star_index[np.where(star_cut)]
  N_stars = np.sum(star_cut)
  print '\tNum. Images Where a star contributes > 600 Times:', N_stars
  print '\tReshape Star Data Matrix from', star_mag_matrix.shape, 'to', star_mag_matrix[star_cut,:].shape
  star_mag_matrix     = star_mag_matrix[star_cut,:]
  star_mag_err_matrix = star_mag_err_matrix[star_cut,:]
  star_mask_matrix    = star_mask_matrix[star_cut,:]
  '''

  ########################################################################################

  # No point in plotting as a function of time.
  '''
  cols = cm.plasma(np.linspace(0, 1, N_stars))
  figure()
  ax = subplot(311)
  for k in range(0,len(star_mags)):
    plot(times, np.ma.masked_where(star_mask_matrix, star_mag_matrix)[k,:], '.', color=cols[k])
  y1,y2 = ax.get_ylim()
  ylim(y2,y1)
  ax1 = subplot(312)
  for k in range(0,len(star_mags)):
    plot(times, time_var + star_mags[k], '.', color=cols[k])
  y1,y2 = ax1.get_ylim()
  ylim(y2,y1)
  ax2 = subplot(313)
  for k in range(0,len(star_mags)):
    plot(times, time_var, 'k.')
  y1,y2 = ax2.get_ylim()
  ylim(y2,y1)
  '''

  ########################################################################################
  # Plot Instrument Magnitudes w/o and w/ time variable corrections.
  cols = cm.plasma(np.linspace(0, 1, N_stars))
  figure(figsize=(10,12))
  sorted_k = np.argsort(star_mags)
  cc=0
  for k in sorted_k: # loop through stars in order of brightness
    cut = star_mask_matrix[k,:]

    ax1=subplot(211)
    errorbar(image_num[cut], star_mag_matrix[k,:][cut] , star_mag_err_matrix[k,:][cut], fmt=',', color=cols[cc])

    ax2=subplot(212)
    errorbar(image_num[cut], star_mag_matrix[k,:][cut] - time_var[cut] , star_mag_err_matrix[k,:][cut], fmt=',', color=cols[cc])
    cc+=1

    med = np.median(star_mag_matrix[k,:][cut] - time_var[cut])
    MAD = 1.4826*np.median(np.abs(star_mag_matrix[k,:][cut] - time_var[cut]-med))
    outlier_cut =  np.abs(star_mag_matrix[k,:] - time_var -med )/MAD > 2.
    print '\t\t%d\t%1.2f\t%1.2f\t%1.2f\t%d'%(k, med,MAD, np.max(np.abs(star_mag_matrix[k,:][cut] - time_var[cut] -med ))/MAD, np.sum(outlier_cut))
    #print 'outlier_cut.shape', outlier_cut.shape
    star_mask_matrix[k,outlier_cut] = False

  ax1=subplot(211)
  y1,y2 = ax1.get_ylim()
  ylim(y2,y1)
  xlabel('Image Index')
  ylabel('Instrument Mag')
  title('Initial Time-Variable Photometric Corrections')
  grid(True)
  ax2=subplot(212)
  grid(True)
  y1,y2 = ax2.get_ylim()
  ylim(y2,y1)
  xlabel('Image Index')
  ylabel('Time-Variability Corrected\nInstrument Mag')
  subplots_adjust(left=0.125)

  ########################################################################################

  print '\n\tRe-Evaluating Initial Estimates After Outlier Cuts'
  star_mags = np.zeros(N_stars)
  for ii in range(0,N_stars):
    star_mags[ii] = np.median(star_mag_matrix[ii,:][star_mask_matrix[ii,:]])
    #star_mags[ii] = np.mean(star_mag_matrix[ii,:][star_mask_matrix[ii,:]])

  time_var = np.zeros(star_mag_matrix.shape[1])
  for jj in range(0, star_mag_matrix.shape[1]):
    res = star_mag_matrix[:,jj] - star_mags[:]
    time_var[jj] = np.median(res[star_mask_matrix[:,jj]])
    #time_var[jj] = np.mean(res[star_mask_matrix[:,jj]])

  print '\n\tCheck Magnitudes and Time Series were all Estimates'
  print '\t\tnp.sum(star_mags==0.)',np.sum(star_mags==0.)
  print '\t\tnp.sum(time_var==0.)',np.sum(time_var==0.)

  ########################################################################################
  # REMOVE TIME POINTS WHERE NO STARS PARTICIPATE

  time_cut = time_var!=0.
  print '\n\tReshape Star Data Matrix from', star_mag_matrix.shape, 'to', star_mag_matrix[:,time_cut].shape

  time_index = np.where(time_cut)[0]
  star_mag_matrix     = star_mag_matrix[:,time_cut]
  star_mag_err_matrix = star_mag_err_matrix[:,time_cut]
  star_mask_matrix    = star_mask_matrix[:,time_cut]
  time_var            = time_var[time_cut]
  times               = times[time_cut]
  image_num           = image_num[time_cut]

  ########################################################################################

  print '\n\tRe-making Initial Estimates of Average Fluxes and Time Variability'
  star_mags = np.zeros(N_stars)
  for ii in range(0,N_stars):
    star_mags[ii] = np.median(star_mag_matrix[ii,:][star_mask_matrix[ii,:]])
    #star_mags[ii] = np.mean(star_mag_matrix[ii,:][star_mask_matrix[ii,:]])

  time_var = np.zeros(star_mag_matrix.shape[1])
  for jj in range(0, star_mag_matrix.shape[1]):
    res = star_mag_matrix[:,jj] - star_mags[:]
    time_var[jj] = np.median(res[star_mask_matrix[:,jj]])
    #time_var[jj] = np.mean(res[star_mask_matrix[:,jj]])

  print '\n\tCheck Magnitudes and Time Series were all Estimates'
  print '\t\tnp.sum(star_mags==0.)',np.sum(star_mags==0.)
  print '\t\tnp.sum(time_var==0.)',np.sum(time_var==0.)


  ########################################################################################

  num_image_contributions_by_star = np.sum(star_mask_matrix, axis=1)
  print 'sorted(num_image_contributions_by_star)',sorted(num_image_contributions_by_star)

  ########################################################################################
  # Plot Instrument Magnitudes w/o and w/ time variable corrections after outlier removal.

  cols = cm.plasma(np.linspace(0, 1, N_stars))
  figure(figsize=(10,12))
  sorted_k = np.argsort(star_mags)
  cc=0
  for k in sorted_k: # loop through stars in order of brightness
    cut = star_mask_matrix[k,:]

    ax1=subplot(211)
    errorbar(image_num[cut], star_mag_matrix[k,:][cut] , star_mag_err_matrix[k,:][cut], fmt=',', color=cols[cc])

    ax2=subplot(212)
    errorbar(image_num[cut], star_mag_matrix[k,:][cut] - time_var[cut] , star_mag_err_matrix[k,:][cut], fmt=',', color=cols[cc])
    cc+=1

    med = np.median(star_mag_matrix[k,:][cut] - time_var[cut])
    MAD = 1.4826*np.median(np.abs(star_mag_matrix[k,:][cut] - time_var[cut]-med))
    #outlier_cut =  np.abs(star_mag_matrix[k,:] - time_var -med )/MAD > 2.
    print '\t\t%d\t%1.2f\t%1.2f\t%1.2f'%(k, med,MAD, np.max(np.abs(star_mag_matrix[k,:][cut] - time_var[cut] -med ))/MAD)
    #print 'outlier_cut.shape', outlier_cut.shape
    #star_mask_matrix[k,outlier_cut] = False

  ax1=subplot(211)
  y1,y2 = ax1.get_ylim()
  ylim(y2,y1)
  xlabel('Image Index')
  ylabel('Instrument Mag')
  title('Initial Time-Variable Photometric Corrections')
  grid(True)
  ax2=subplot(212)
  grid(True)
  y1,y2 = ax2.get_ylim()
  ylim(y2,y1)
  xlabel('Image Index')
  ylabel('Time-Variability Corrected\nInstrument Mag')
  subplots_adjust(left=0.15)


  ######################################################################################## 
  '''

  print '\n\tCompare results to SDSS Magnitudes '
  global_ZP = np.median(np.array(star_table['SDSS_r'])[star_index]-star_mags[:N_stars])
  global_ZP_MAD = 1.4826*np.median(np.abs(np.array(star_table['SDSS_r'])[star_index]-star_mags[:N_stars] - global_ZP)) 
  print '\t\tEstimated global_ZP %1.2f +/- %1.2f'%(global_ZP, global_ZP_MAD)
  sorted_k = np.argsort(star_mags[:N_stars])
  print '\t\tStar\tm_I\tm_SDSS\tDiff\tglobal_ZP Res'
  for k in sorted_k:
    print '\t\t%d\t%1.2f\t%1.2f\t%1.2f\t%+1.2f'%( k, star_mags[k], np.array(star_table['SDSS_r'])[star_index[k]], np.array(star_table['SDSS_r'])[star_index[k]] - star_mags[k], np.array(star_table['SDSS_r'])[star_index[k]] - star_mags[k]-global_ZP)

  ######################################################################################## 
  print '\n\tCut Star With Deviations Greater than 3 times the SDSS Flux Error'
  ZP_dev = np.array(star_table['SDSS_r'])[star_index] - star_mags - global_ZP
  flux_err = np.array(star_table['FLUXERR_AUTO_1'])[star_index]/np.array(star_table['FLUX_AUTO_1'])[star_index]
  cal_cut = np.abs( ZP_dev / flux_err ) < 5.0
  #print 'cal_cut', cal_cut
  star_index = star_index[np.where(cal_cut)]
  #print 'star_index', star_index
  print '\t\tnp.sum(cal_cut)', np.sum(cal_cut)
  
  N_stars = np.sum(cal_cut)
  #print '\t\tstar_mag_matrix.shape, cal_cut.shape ',star_mag_matrix.shape, cal_cut.shape 
  print '\t\tstar_mag_matrix[cal_cut,:].shape, cal_cut.shape', star_mag_matrix[cal_cut,:].shape, cal_cut.shape 
  star_mag_matrix     = star_mag_matrix[cal_cut,:]
  star_mag_err_matrix = star_mag_err_matrix[cal_cut,:]
  star_mask_matrix    = star_mask_matrix[cal_cut,:]
  star_mags           = star_mags[cal_cut]

  ######################################################################################## 

  print '\n\tCompare results to SDSS Magnitudes '
  global_ZP = np.median(np.array(star_table['SDSS_r'])[star_index]-star_mags[:N_stars])
  global_ZP_MAD = 1.4826*np.median(np.abs(np.array(star_table['SDSS_r'])[star_index]-star_mags[:N_stars] - global_ZP)) 
  print '\t\tEstimated global_ZP %1.2f +/- %1.2f'%(global_ZP, global_ZP_MAD)
  sorted_k = np.argsort(star_mags[:N_stars])
  print '\t\tStar\tm_I\tm_SDSS\tDiff\tglobal_ZP Res'
  for k in sorted_k:
    print '\t\t%d\t%1.2f\t%1.2f\t%1.2f\t%+1.2f'%( k, star_mags[k], np.array(star_table['SDSS_r'])[star_index[k]], np.array(star_table['SDSS_r'])[star_index[k]] - star_mags[k], np.array(star_table['SDSS_r'])[star_index[k]] - star_mags[k]-global_ZP)
  '''
  ######################################################################################## 
  print '\n\tBegin Fit'

  from scipy.optimize import curve_fit
  print '\t\tlen(star_mags), len(time_var)', len(star_mags), len(time_var)
  p0 = np.concatenate([star_mags, time_var])
  print '\t\tp0.shape', p0.shape

  ref_mag = p0[0]
  p0 -= ref_mag
  #print len(p0), len(star_mag_matrix[star_mask_matrix])

  # create a function to fit to 
  # this is meant to fit an array of time series to a set of means plus one time-varying component
  def func((i2d,t2d), *parms):
      N = len(t2d)
      M = len(i2d)
      star_mags = np.array(parms[0:np.max(i2d)+1])
      time_var = np.array(parms[np.max(i2d)+1:])
      _mags = np.zeros(N)
      for k in range(0,N):
         _mags[k] = star_mags[i2d[k]] + time_var[t2d[k]]
      return _mags.ravel()

  
  def func_fixed_ref((i2d,t2d), *parms):
      # insert line so that the first parm is always 0.
      parms = np.concatenate([[0.], parms])
      #print parms
      return func((i2d,t2d), *parms)

  print 'START CURVE_FIT'

  I2d, T2d = np.indices(star_mag_matrix.shape)
  popt, pcov = curve_fit(func_fixed_ref, 
                        (I2d[star_mask_matrix], T2d[star_mask_matrix]), 
                        star_mag_matrix[star_mask_matrix], 
                        sigma=star_mag_err_matrix[star_mask_matrix], 
                        p0=p0[1:], # reference is excluded
                        check_finite=True, 
                        #method = 'lm', # only works with unbounded problems
                        method = 'trf',
                        absolute_sigma=True)

  #exit()
  popt = np.concatenate([[0.], popt])
  perr = np.concatenate([[0.],np.sqrt(np.diag(pcov))])

  popt[:N_stars] += ref_mag
  popt[N_stars:] -= ref_mag
  #print '\t\tpopt', popt
  #print '\t\tnp.diag(pcov)', perr

  print '\tpopt.shape',popt.shape 
  print '\tpcov.shape',perr.shape 

  fit_mag          = popt[:N_stars]
  fit_time_var     = popt[N_stars:]
  fit_mag_err      = perr[:N_stars]
  fit_time_var_err = perr[N_stars:]

  fnm_list = np.array(fnm_list)
  np.savez('phot_cal.npz', time_index = time_index,
                           star_index = star_index,
                           image_num  = image_num,
                           fnm_list  = fnm_list[time_index],
                           fit_time_var = fit_time_var,
                           fit_time_var_err = fit_time_var_err,
                           fit_mag = fit_mag,
                           fit_mag_err = fit_mag_err,
                           #global_ZP = global_ZP,
                           star_mag_matrix = star_mag_matrix,
                           star_mag_err_matrix = star_mag_err_matrix,
                           star_mask_matrix = star_mask_matrix
                           )
  ######################################################################################## 

  figure(figsize=(10,12))
  ax=subplot(211)
  errorbar(image_num, fit_time_var, yerr=fit_time_var_err, fmt='.')
  ylabel('Time-Variable Magnitude Correction, mag')
  grid(True)
  y1,y2 = ax.get_ylim()
  ylim(y2,y1)
  title('Best Fit Time Variability')
  subplot(212)
  semilogy(image_num, fit_time_var_err, '.')
  xlabel('Image Index')
  grid(True, which='both')
  ylabel('Uncertainty, mag')
  print '\ttimes.shape',times.shape
  print '\tstar_mag_matrix.shape',star_mag_matrix.shape
  print '\tstar_mag_err_matrix.shape',star_mag_err_matrix.shape

  figure(figsize=(10,8))
  ax=subplot(111)
  sorted_k = np.argsort(fit_mag)
  cc=0
  for k in sorted_k:
    cut = star_mask_matrix[k,:]
    errorbar(image_num[cut], star_mag_matrix[k,:][cut] - fit_time_var[cut] , star_mag_err_matrix[k,:][cut], fmt='.', color=cols[cc])
    cc+=1
  y1,y2 = ax.get_ylim()
  ylim(y2,y1)
  xlabel('Image Index')
  ylabel('Time Variability Corrected\nInstrument Magnitude')
  subplots_adjust(left=0.15)
  title('Fitted Variability Corrections')
  grid(True)

  '''
  print '\n\tCompare results to SDSS Magnitudes '
  global_ZP = np.median(np.array(star_table['SDSS_r'])[star_index]-fit_mag[:N_stars])
  global_ZP_MAD = 1.4826*np.median(np.abs(np.array(star_table['SDSS_r'])[star_index]-fit_mag[:N_stars] - global_ZP)) 
  print '\t\tEstimated global_ZP %1.2f +/- %1.2f'%(global_ZP, global_ZP_MAD)
  sorted_k = np.argsort(fit_mag[:N_stars])
  print '\t\tStar\tm_I\tm_SDSS\tDiff\tglobal_ZP Res\terr'
  for k in sorted_k:
    cut = star_mask_matrix[k,:]
    star_ZP = np.array(star_table['SDSS_r'])[star_index[k]] - fit_mag[k]
    res  =  star_ZP-global_ZP
    bias = np.median(star_mag_matrix[k,:][cut] - fit_time_var[cut] + global_ZP - np.array(star_table['SDSS_r'][star_index[k]]))
    err = 1.4826*np.median(np.abs(star_mag_matrix[k,:][cut] - fit_time_var[cut] + global_ZP - np.array(star_table['SDSS_r'][star_index[k]]) - bias))
    print '\t\t%d\t%1.2f\t%1.2f\t%1.2f\t%+1.3f\t\t%+1.3f\t%+1.3f'%( k, fit_mag[k], np.array(star_table['SDSS_r'])[star_index[k]], star_ZP, res, err, np.array(star_table['SDSS_g-r'])[star_index[k]])
    #print '\t\t%d\t%1.2f\t%1.2f\t%1.2f\t%+1.4f\t%1.5f\t%1.5f'%( k, fit_mag[k], np.array(star_table['SDSS_r'])[star_index[k]], star_ZP, res,  bias, err)


  figure()
  ax=subplot(111)
  #ax.set_yscale('log')
  sorted_k = np.argsort(fit_mag)
  #for k in range(0,N_stars):
  cc=0
  #for k in range(0, len(fit_mag)):
  for k in sorted_k:
    #cut = star_mag_matrix[k,:]>0.
    cut = star_mask_matrix[k,:]
    #print k, times.shape, star_mag_matrix[k,:].shape 
    #plot(times[cut], star_mag_matrix[k,:][cut], '.-')
    #errorbar(times[cut], star_mag_matrix[k,:][cut], star_mag_err_matrix[k,:][cut], fmt='.')
    #errorbar(times[cut], star_mag_matrix[k,:][cut] - fit_time_var[cut] , star_mag_err_matrix[k,:][cut], fmt='.')
    errorbar(image_num[cut], star_mag_matrix[k,:][cut] - fit_time_var[cut] + global_ZP, star_mag_err_matrix[k,:][cut], fmt='.', color=cols[cc])
    plot(image_num[cut], np.array(star_table['SDSS_r'])[star_index[k]]*np.ones(len(image_num[cut])), '--', color = cols[cc])
    cc+=1
  y1,y2 = ax.get_ylim()
  ylim(y2,y1)
  ylabel('Calibrated Fluxes, mag')
  xlabel('Image Index')
  '''
  return popt[N_stars:], perr[N_stars:]
#Star_chi_max_cut = 4.5
#Star_chi_sq_cut = 1.15
#num_Star_good_cut = 3
#Star_fractional_uncertainty_cut = 0.3

######################################################################################################
#['nmax', 'shapelet_coeffs', 'qsr_wg_parms', 'star_S_CCD_unc', 'star_index_list', 'inputFile', 'mmax', 'APASS_sig_S_CCD', 'APASS_max_chi', 'outFileTag', 'beta_mean', 'star_max_chi', 'qsr_wg_covar', 'filter', 'qsr_chisq', 'star_S_CCD', 'APASS_index_list', 'APASS_S_CCD', 'mjd_obs', 'star_chi_sq', 'APASS_chi_sq', 'beta_unc', 'qsr_max_chi', 'readnoise']

def flux_corrected_qsr_curves(fnm_list, rel_phot, rel_phot_err, kind='wg', t0=57008.):
  # index the times, magnitudes, errors, and create mask matrix for the values that pass our cuts.
  times   = [] 
  img_1   = []
  img_2   = []
  img_3   = []
  img_4   = []
  img_gal = []
  img_1_err = []
  img_2_err = []
  img_3_err = []
  img_4_err = []
  img_gal_err = []

  for i in range(0,len(fnm_list)):
    f = np.load(fnm_list[i])
    qual_cut = np.all(np.sqrt(np.diagonal(f['qsr_%s_covar'%kind]))[10:14] / f['qsr_%s_parms'%kind][10:14]<QSR_fractional_uncertainty_cut)

    qual_cut = np.logical_and(f['qsr_%s_chisq'%kind]<QSR_chi_sq_cut, f['qsr_%s_max_chi'%kind]<QSR_chi_max_cut)
    qual_cut = np.logical_and(np.all(np.sqrt(np.diagonal(f['qsr_%s_covar'%kind]))[10:14] / f['qsr_%s_parms'%kind][10:14]<QSR_fractional_uncertainty_cut), qual_cut)
    qual_cut = np.logical_and(np.all(f['qsr_%s_parms'%kind][10:14]>0.), qual_cut)
    qual_cut = np.logical_and(np.all(np.sqrt(np.diagonal(f['qsr_%s_covar'%kind]))[10:14]>0.), qual_cut)
    qual_cut = np.logical_and(np.all(np.isfinite(f['qsr_%s_parms'%kind][10:14])), qual_cut)
    qual_cut = np.logical_and(np.all(np.isfinite(np.sqrt(np.diagonal(f['qsr_%s_covar'%kind]))[10:14])), qual_cut)
    if (not qual_cut):
        print qual_cut
    if(qual_cut):
        times.append(f['mjd_obs']-t0)
        img_1.append(AL.flux2magnitude(f['qsr_%s_parms'%kind][10]))
        img_2.append(AL.flux2magnitude(f['qsr_%s_parms'%kind][11]))
        img_3.append(AL.flux2magnitude(f['qsr_%s_parms'%kind][12]))
        img_4.append(AL.flux2magnitude(f['qsr_%s_parms'%kind][13]))
        #img_gal.append(AL.flux2magnitude(f['qsr_%s_parms'%kind][17]))
        img_1_err.append(AL.fluxErr2magErr(f['qsr_%s_parms'%kind][10],np.sqrt(np.diagonal(f['qsr_%s_covar'%kind]))[10]))
        img_2_err.append(AL.fluxErr2magErr(f['qsr_%s_parms'%kind][11],np.sqrt(np.diagonal(f['qsr_%s_covar'%kind]))[11]))
        img_3_err.append(AL.fluxErr2magErr(f['qsr_%s_parms'%kind][12],np.sqrt(np.diagonal(f['qsr_%s_covar'%kind]))[12]))
        img_4_err.append(AL.fluxErr2magErr(f['qsr_%s_parms'%kind][13],np.sqrt(np.diagonal(f['qsr_%s_covar'%kind]))[13]))
        #img_gal_err.append(AL.fluxErr2magErr(f['qsr_%s_parms'%kind][17], np.sqrt(np.diagonal(f['qsr_%s_covar'%kind]))[17]))
    f.close()
  times = np.array(times)
  img_1 = np.array(img_1)
  img_2 = np.array(img_2)
  img_3 = np.array(img_3)
  img_4 = np.array(img_4)
  img_gal = np.array(img_gal)
  img_1_err = np.array(img_1_err)
  img_2_err = np.array(img_2_err)
  img_3_err = np.array(img_3_err)
  img_4_err = np.array(img_4_err)
  img_gal_err = np.array(img_gal_err)


  figure()
  ax1 = subplot(411)
  errorbar(times, img_1, yerr=img_1_err, fmt = '.')
  grid(True)
  y11, y12 = ax1.get_ylim()
  ylim(y12,y11)
  ax2 = subplot(412, sharex=ax1, sharey=ax1)
  errorbar(times, img_2, yerr=img_2_err, fmt = '.')
  grid(True)
  y21, y22 = ax2.get_ylim()
  ylim(y22,y21)
  ax3 = subplot(413, sharex=ax1, sharey=ax1)
  errorbar(times, img_3, yerr=img_3_err, fmt = '.')
  grid(True)
  y31, y32 = ax3.get_ylim()
  ylim(y32,y31)
  ax4 = subplot(414, sharex=ax1, sharey=ax1)
  errorbar(times, img_4, yerr=img_4_err, fmt = '.')
  grid(True)
  y41, y42 = ax4.get_ylim()
  ylim(y42,y41)
  #ax5 = subplot(515)
  #errorbar(times, img_gal, yerr=img_gal_err, fmt = '.')
  #grid(True)
  #y1, y2 = ax5.get_ylim()
  #ylim(y2,y1)

  img_1_rel_mean = np.average(img_1-rel_phot, weights = 1./(img_1_err**2+rel_phot_err**2))
  img_2_rel_mean = np.average(img_2-rel_phot, weights = 1./(img_2_err**2+rel_phot_err**2))
  img_3_rel_mean = np.average(img_3-rel_phot, weights = 1./(img_3_err**2+rel_phot_err**2))
  img_4_rel_mean = np.average(img_4-rel_phot, weights = 1./(img_4_err**2+rel_phot_err**2))

  figure()
  ax1 = subplot(411)
  errorbar(times, img_1-rel_phot - img_1_rel_mean, yerr=np.sqrt(img_1_err**2+rel_phot_err**2), fmt = '.')
  grid(True)
  y11, y12 = ax1.get_ylim()
  ylim(y12,y11)
  ax2 = subplot(412, sharex=ax1, sharey=ax1)
  errorbar(times, img_2-rel_phot - img_2_rel_mean, yerr=np.sqrt(img_2_err**2+rel_phot_err**2), fmt = '.')
  grid(True)
  y21, y22 = ax2.get_ylim()
  ylim(y22,y21)
  ax3 = subplot(413, sharex=ax1, sharey=ax1)
  errorbar(times, img_3-rel_phot - img_3_rel_mean, yerr=np.sqrt(img_3_err**2+rel_phot_err**2), fmt = '.')
  grid(True)
  y31, y32 = ax3.get_ylim()
  ylim(y32,y31)
  ax4 = subplot(414, sharex=ax1, sharey=ax1)
  errorbar(times, img_4-rel_phot - img_4_rel_mean, yerr=np.sqrt(img_4_err**2+rel_phot_err**2), fmt = '.')
  grid(True)
  y41, y42 = ax4.get_ylim()
  ylim(y42,y41)
  #ax5 = subplot(515)
  #errorbar(times, img_gal-rel_phot, yerr=np.sqrt(img_gal_err**2+rel_phot_err**2), fmt = '.')
  #grid(True)
  #y1, y2 = ax5.get_ylim()
  #ylim(y2,y1)

  bins = np.arange(0., 400. + 3./24., 3./24.)
  figure()
  ax1 = subplot(411)
  binned_t, m1, e1 = binned_light_curve(times, img_1 - rel_phot - img_1_rel_mean, np.sqrt(img_1_err**2 + rel_phot_err**2), bins)
  errorbar(binned_t, m1, yerr=e1, fmt='b.')
  #plt.gca().invert_yaxis()
  ylim(y12,y11)
  grid(True)

  ax2 = subplot(412, sharex=ax1, sharey=ax1)
  binned_t, m2, e2 = binned_light_curve(times, img_2 - rel_phot - img_2_rel_mean, np.sqrt(img_2_err**2 + rel_phot_err**2), bins)
  errorbar(binned_t, m2, yerr=e2, fmt='b.')
  #plt.gca().invert_yaxis()
  ylim(y22,y21)
  grid(True)

  ax3 = subplot(413, sharex=ax1, sharey=ax1)
  binned_t, m3, e3 = binned_light_curve(times, img_3 - rel_phot - img_3_rel_mean, np.sqrt(img_3_err**2 + rel_phot_err**2), bins)
  errorbar(binned_t, m3, yerr=e3, fmt='b.')
  #plt.gca().invert_yaxis()
  ylim(y32,y31)
  grid(True)

  ax4 = subplot(414, sharex=ax1, sharey=ax1)
  binned_t, m4, e4 = binned_light_curve(times, img_4 - rel_phot - img_4_rel_mean, np.sqrt(img_4_err**2 + rel_phot_err**2), bins)
  errorbar(binned_t, m4, yerr=e4, fmt='b.')
  #plt.gca().invert_yaxis()
  ylim(y42,y41)
  grid(True)

  return times, img_1 - rel_phot, img_2 - rel_phot, img_3 - rel_phot, img_4 - rel_phot, np.sqrt(img_1_err**2 + rel_phot_err**2), np.sqrt(img_2_err**2 + rel_phot_err**2), np.sqrt(img_3_err**2 + rel_phot_err**2), np.sqrt(img_4_err**2 + rel_phot_err**2)
    
######################################################################################################


fnm_list = []
number_found = 0
number_with_warning = 0
number_not_found = 0

import glob
psf_dir = '/halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/Phot'
#files = glob.glob(psf_dir+'/*.npz')
#files = sorted(files)

fnm_list = []
for k in range(0,1458):
  fnm = psf_dir+'/image_%04d_Phot.npz'%k # first and second run combined
  if(os.path.isfile(fnm)):
    number_found += 1
    fnm_list.append(fnm)
    f = np.load(fnm)
    #print k, f.keys()
    f.close()
  #'''
  if(not os.path.isfile(fnm)):
    number_not_found += 1
  #'''
print ''
print 'number found\t\t', number_found
#print '\tnumber with a warning', number_with_warning
print 'number not found\t', number_not_found
print 'total\t\t\t', number_found + number_not_found
print ''
t0 = 56972.

#plot_APASS_chi(fnm_list)
print 'number_found\t', len(fnm_list)

fnm_list_f = filter_cuts(fnm_list, filter = 'rp')
print 'number red:\t', len(fnm_list_f)
fnm_list_1 = PSF_quality_cuts(fnm_list_f)
print 'good PSF\t', len(fnm_list_1)

fnm_list_2 = Star_quality_cuts(fnm_list_1)
print 'good Star\t', len(fnm_list_2)
#fnm_list_3 = QSR_quality_cuts(fnm_list_2, kind = qsr_fit_kind)
#print 'good QSR\t', len(fnm_list_3)

rel_phot, rel_phot_err = flux_corrected_curves(fnm_list_2, t0=t0)
show()
exit()

t, m1, m2, m3, m4, e1, e2, e3, e4 = flux_corrected_qsr_curves(fnm_list_3, rel_phot, rel_phot_err, kind = qsr_fit_kind, t0=t0)

np.savez('qsr_lightcurve_%s.npz'%qsr_fit_kind, time=t, mag_1=m1, mag_2=m2, mag_3=m3, mag_4=m4, err_1=e1, err_2=e2, err_3=e3, err_4=e4 )
# SHOW DISTRIBUTIONS FOR ALL SURVIVING IMAGES
#APASS_quality_cuts(fnm_list_3)
#Star_quality_cuts(fnm_list_3)
#QSR_quality_cuts(fnm_list_3)

#print 'FAILS FIT'
#print set(fnm_list) - set(fnm_list_1)
show()


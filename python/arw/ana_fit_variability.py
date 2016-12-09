
from pylab import *
import os
import AnalysisLibrary as AL
rcParams['figure.facecolor'] = 'white'
rcParams['font.size'] = 16
from astropy.io import ascii
'''
np.savez('phot_cal.npz', time_index = time_index,
                       star_index = star_index,
                       fnm_list  = fnm_list[time_index],
                       fit_time_var = fit_time_var,
                       fit_time_var_err = fit_time_var_err,
                       fit_mag = fit_mag,
                       fit_mag_err = fit_mag_err,
                       global_ZP = global_ZP
                       )
'''

global_zero_point = 29.425
chi_sq_cut = 1.3
max_chi_cut = 4.3
frac_unc_M1_M2 = 0.3

def absolute_cal(npz_filename):
    f = np.load(npz_filename)
    #print f.keys()
    #print ''
    image_num = f['image_num']
    star_table = ascii.read('sample.csv')
    fit_mag = f['fit_mag']
    star_index = f['star_index']
    star_mag_matrix = f['star_mag_matrix']
    star_mag_err_matrix = f['star_mag_err_matrix']
    star_mask_matrix = f['star_mask_matrix']
    fit_time_var = f['fit_time_var']
    f.close()
    print 'len(fit_mag)', len(fit_mag)
    N_stars = len(fit_mag)

    SDSS_r = np.array(star_table['SDSS_r'])[star_index]

    print '\n\tCompare results to SDSS Magnitudes '
    flux_err = np.array(star_table['FLUXERR_AUTO_1'])[star_index] / np.array(star_table['FLUX_AUTO_1'])[star_index]
    star_ZP = SDSS_r - fit_mag
    global_ZP = np.median(star_ZP)
    ZP_res = star_ZP - global_ZP
    global_ZP_MAD = 1.4826*np.median(np.abs(ZP_res))
    print '\t\tEstimated global_ZP %1.3f +/- %1.3f'%(global_ZP, global_ZP_MAD)

    ZP_res_cut = np.abs((ZP_res) / flux_err)<3.

    global_ZP = np.median(star_ZP[ZP_res_cut])
    global_ZP_MAD = 1.4826*np.median(np.abs(ZP_res[ZP_res_cut]))

    ZP_res = star_ZP - global_ZP
    print '\t\tEstimated global_ZP %1.3f +/- %1.3f'%(global_ZP, global_ZP_MAD)

    figure(figsize=(10,8))
    sorted_k = np.argsort(fit_mag)
    print '\t\tStar\tm_I\tm_SDSS\tDiff\tglobal_ZP_Res\terr\tcolor\tDev_sig'
    for k in sorted_k:
      if(ZP_res_cut[k]):
        cut = star_mask_matrix[k,:]
        bias = np.median(star_mag_matrix[k,:][cut] - fit_time_var[cut] + global_ZP - SDSS_r[k])
        err = 1.4826*np.median(np.abs(star_mag_matrix[k,:][cut] - fit_time_var[cut] + global_ZP - SDSS_r[k] - bias))
        print '\t\t%d\t%1.2f\t%1.2f\t%1.2f\t%+1.3f\t\t%+1.3f\t%+1.3f\t%+1.3f'%( k, fit_mag[k], SDSS_r[k], star_ZP[k], ZP_res[k], err, np.array(star_table['SDSS_g-r'])[star_index[k]], ZP_res[k]/flux_err[k])

        subplot(221)
        plot([np.array(star_table['SDSS_g-r'])[star_index[k]]],[ZP_res[k]], 'ko')
        xlabel('SDSS g-r, mag')
        ylabel('ZP deviation, mag')
        ylim(-0.2, 0.2)
        grid(True)
        subplot(222)
        plot([np.array(star_table['SDSS_g-r'])[star_index[k]]],[err], 'ko')
        xlabel('SDSS g-r, mag')
        ylabel('Std. Dev., mag')
        grid(True)
        ylim(0., 0.05)
        subplot(223)
        plot([flux_err[k]],[ZP_res[k]], 'ko')
        xlabel('SDSS_r error, mag')
        ylabel('ZP deviation, mag')
        grid(True)
        ylim(-0.2, 0.2)
        xlim(0.,0.2)
       # xticks(rotation=45)
        subplot(224)
        plot([flux_err[k]],[err], 'ko')
        xlabel('SDSS_r error, mag')
        ylabel('Std. Dev., mag')
        grid(True)
        ylim(0., 0.05)
        xlim(0.,0.2)
        #xticks(rotation=45)
    subplots_adjust(hspace=0.25, wspace=0.4, bottom=0.1)


    '''
    figure(figsize=(10,8))
    subplot(221)
    plot( np.array(flux_err), np.array(ZP_res) , 'ko' )
    subplot(222)
    plot( np.array(flux_err), np.array(ZP_res)/np.array(flux_err) , 'ko' )
    subplot(223)
    plot( np.array(ZP_res), np.array(ZP_res)/np.array(flux_err) , 'ko' )
    '''

    figure(figsize=(10,8))
    subplot(211)
    hist(ZP_res, bins=np.linspace(-0.2,0.2,41))
    xlim(-0.2, 0.2)
    xlabel('Star ZP - Median ZP, mag')
    grid(True)
    subplot(212)
    hist(np.array(ZP_res) / np.array(flux_err), bins=np.linspace(-10.,10.,41))
    xlim(-10.,10.)
    xticks(np.linspace(-10.,10.,11) )
    xlabel(' (Star ZP - Median ZP) / SDSS Flux Error')
    grid(True)


    print np.sort(np.abs(np.array(ZP_res) / np.array(flux_err)))[::-1]
    #cols = cm.plasma(np.linspace(0, 1, N_stars))
    cols = cm.plasma(np.linspace(0, 1, np.sum(ZP_res_cut)))


    figure()
    ax=subplot(111)
    #ax.set_yscale('log')
    sorted_k = np.argsort(fit_mag)
    #for k in range(0,N_stars):
    cc=0
    #for k in range(0, len(fit_mag)):
    print global_ZP
    for k in sorted_k:
        if(not ZP_res_cut[k]): continue
        #cut = star_mag_matrix[k,:]>0.
        cut = star_mask_matrix[k,:]
        #errorbar(image_num[cut], star_mag_matrix[k,:][cut] - fit_time_var[cut] + global_ZP, star_mag_err_matrix[k,:][cut], fmt='.', color=cols[cc])
        errorbar(image_num[cut], star_mag_matrix[k,:][cut]  - fit_time_var[cut]  + global_ZP, star_mag_err_matrix[k,:][cut], fmt='.', color=cols[cc])
        #plot(image_num[cut], np.array(star_table['SDSS_r'])[star_index[k]]*np.ones(len(image_num[cut])), '--', color = cols[cc])
        #errorbar(image_num[cut], star_mag_matrix[k,:][cut] - fit_time_var[cut] + global_ZP, star_mag_err_matrix[k,:][cut], fmt='.', color=cols[cc])
        plot([0.,1600.], [SDSS_r[k], SDSS_r[k]], '--', color = cols[cc])
        cc+=1
    y1,y2 = ax.get_ylim()
    ylim(y2,y1)
    ylabel('Calibrated Fluxes, mag')
    xlabel('Image Index')

######################################################

def qsr_file_quality_parameters(qsr_fnm):
  g = np.load(qsr_fnm)
  #print g.keys()
  #print g['parm_list']
  chi_sq, max_chi =  g['chi_sq'], g['max_chi']
  popt_wg = g['popt_wg']
  pcov_wg = g['pcov_wg']
  #x0,y0, S_CCD_0, S_CCD_1, S_CCD_2, S_CCD_3, N_bkg, S_CCD_LensGal = g['popt_wg']
  #x0_err,y0_err, S_CCD_0_err, S_CCD_1_err, S_CCD_2_err, S_CCD_3_err, N_bkg, S_CCD_LensGal = g['popt_wg']
  g.close()
  return chi_sq, max_chi, np.array(popt_wg), np.array(np.sqrt(np.diag(pcov_wg))) 

######################################################

def plot_qsr_chi(chi_sq_list, max_chi_list):

  display_cut = np.isfinite(chi_sq_list)
  #display_cut = np.logical_and( np.isfinite(frac_uncertainty), display_cut )
  #display_cut = np.logical_and( frac_uncertainty>0., display_cut )

  figure(figsize=(10,8))
  suptitle('Quad Fit $\chi^2$ and max$|\chi|$')
  h,b = np.histogram( np.log10(chi_sq_list[display_cut]), bins=500)
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
  plot([chi_sq_cut, chi_sq_cut],[0.5, 2.*np.max(h)], 'r--', lw=2)
  grid(True, which='both')

  h,b = np.histogram( np.log10(max_chi_list[display_cut]), bins=500)
  ax2 = subplot(224)
  ax2.set_yscale('log')
  ax2.set_xscale('log')
  plot(10**b[1:], h+1.e-3, drawstyle='steps')
  ylim(0.5, 2.*np.max(h))
  plot([max_chi_cut, max_chi_cut],[0.5, 2.*np.max(h)], 'r--', lw=2)
  x21, x22 = ax2.get_xlim()
  x21 = 2.
  xlim(x21, x22)
  grid(True, which='both')
  xlabel(r'max$|\chi|$')

  ax3 = subplot(223, sharex=ax1)
  nbins = 100
  H, xedges, yedges = np.histogram2d(np.log10(chi_sq_list[display_cut]), np.log10(max_chi_list[display_cut]), bins=nbins)
  H = np.rot90(H)
  H = np.flipud(H) 
  Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
  from matplotlib.colors import LogNorm
  ax3.pcolormesh(10**xedges,10**yedges,Hmasked, cmap='viridis', norm=LogNorm()) 
  xlabel(r'$\chi^2$')
  ylabel(r'max$|\chi|$')
  xlim(x11,x12)
  ylim(x21,x22)
  plot([chi_sq_cut, chi_sq_cut], [x21, x22], 'r--', lw=2)
  plot([x11, x12], [max_chi_cut, max_chi_cut], 'r--', lw=2)
  ax3.set_yscale('log')
  ax3.set_xscale('log')


######################################################
def plot_parms(qsr_parm,qsr_parm_err, galFit=True):
    print qsr_parm.shape
    print qsr_parm_err.shape
    figure()
    errorbar(qsr_parm[:,0], qsr_parm[:,1], xerr=qsr_parm_err[:,0], yerr=qsr_parm_err[:,0], fmt='.' )
    xlim(-15.,15.)
    ylim(-15.,15.)
    xlabel('QSR Image x0')
    ylabel('QSR Image y0')

    if(galFit==True):
        print 'Greater than zero amps        ', np.sum(qsr_parm[:,2]>0.), np.sum(qsr_parm[:,3]>0.), np.sum(qsr_parm[:,4]>0.), np.sum(qsr_parm[:,5]>0.), np.sum(qsr_parm[:,6]>0.), np.sum(qsr_parm[:,7]>0.)
        print 'Finite amps                   ', np.sum(np.isfinite(qsr_parm[:,2])), np.sum(np.isfinite(qsr_parm[:,3])), np.sum(np.isfinite(qsr_parm[:,4])), np.sum(np.isfinite(qsr_parm[:,5])), np.sum(np.isfinite(qsr_parm[:,6])), np.sum(np.isfinite(qsr_parm[:,7]))

        print 'Greater than zero amps errors ', np.sum(qsr_parm_err[:,2]>0.), np.sum(qsr_parm_err[:,3]>0.), np.sum(qsr_parm_err[:,4]>0.), np.sum(qsr_parm_err[:,5]>0.), np.sum(qsr_parm_err[:,6]>0.), np.sum(qsr_parm_err[:,7]>0.)
        print 'Finite amps errors            ', np.sum(np.isfinite(qsr_parm_err[:,2])), np.sum(np.isfinite(qsr_parm_err[:,3])), np.sum(np.isfinite(qsr_parm_err[:,4])), np.sum(np.isfinite(qsr_parm_err[:,5])), np.sum(np.isfinite(qsr_parm_err[:,6])), np.sum(np.isfinite(qsr_parm_err[:,7]))

    if(galFit==False):
        print 'Greater than zero amps        ', np.sum(qsr_parm[:,2]>0.), np.sum(qsr_parm[:,3]>0.), np.sum(qsr_parm[:,4]>0.), np.sum(qsr_parm[:,5]>0.), np.sum(qsr_parm[:,6]>0.)
        print 'Finite amps                   ', np.sum(np.isfinite(qsr_parm[:,2])), np.sum(np.isfinite(qsr_parm[:,3])), np.sum(np.isfinite(qsr_parm[:,4])), np.sum(np.isfinite(qsr_parm[:,5])), np.sum(np.isfinite(qsr_parm[:,6]))

        print 'Greater than zero amps errors ', np.sum(qsr_parm_err[:,2]>0.), np.sum(qsr_parm_err[:,3]>0.), np.sum(qsr_parm_err[:,4]>0.), np.sum(qsr_parm_err[:,5]>0.), np.sum(qsr_parm_err[:,6]>0.)
        print 'Finite amps errors            ', np.sum(np.isfinite(qsr_parm_err[:,2])), np.sum(np.isfinite(qsr_parm_err[:,3])), np.sum(np.isfinite(qsr_parm_err[:,4])), np.sum(np.isfinite(qsr_parm_err[:,5])), np.sum(np.isfinite(qsr_parm_err[:,6]))

    figure()
    subplot(211)
    plot(qsr_parm[:,2], '.')
    plot(qsr_parm[:,3], '.')
    plot(qsr_parm[:,4], '.')
    plot(qsr_parm[:,5], '.')
    subplot(212)
    plot(qsr_parm_err[:,2], '.')
    plot(qsr_parm_err[:,3], '.')
    plot(qsr_parm_err[:,4], '.')
    plot(qsr_parm_err[:,5], '.')

    if(galFit):
        figure()
        subplot(211)
        plot(qsr_parm[:,7], '.')
        subplot(212)
        plot(qsr_parm_err[:,7], '.')

    figure()
    subplot(211)
    plot(qsr_parm[:,6], '.')
    subplot(212)
    plot(qsr_parm_err[:,6], '.')

    figure(figsize=(8,12))
    suptitle('Fractional Uncertainties', fontsize=20)
    cc=0
    lab = ['M1', 'S1', 'M2', 'S2', 'Galaxy']
    for k in [2,3,4,5,7]:
        if(not galFit and k==7): 
            continue
        h,b = np.histogram( np.log10(qsr_parm_err[:,k]/qsr_parm[:,k]), bins=linspace(-4,3))
        ax1 = subplot(5,1,cc+1)
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        plot(10**b[1:], h+1.e-3, drawstyle='steps', lw=2, label=lab[cc])
        y11 = 0.5
        y12 = 2.*np.max(h)
        ylim(y11, y12)
        #if (k==2 or k==4):
        plot([frac_unc_M1_M2,frac_unc_M1_M2], [y11,y12], 'r--', lw=2)
        cc+=1
        xlim(1.e-3, 1.e3)
        legend(loc=1)
        grid(True, which='both')
    subplots_adjust(top=0.95)
    xlabel('Fractional Uncertainty on Flux')
    #exit()

def mag_err_and_cut(_flux, _flux_err, _qsr_image_num, _fit_time_var, _fit_time_var_err):
    cut = _flux_err/_flux<0.3
    flux = _flux[cut]
    flux_err = _flux_err[cut]
    qsr_image_num_cut = _qsr_image_num[cut]
    fit_time_var_cut  = _fit_time_var[cut]
    fit_time_var_err_cut  = _fit_time_var_err[cut]
    mag = AL.flux2magnitude(flux)
    mag_err = AL.fluxErr2magErr(flux, flux_err)
    return mag, mag_err, flux, flux_err, qsr_image_num_cut, fit_time_var_cut, fit_time_var_err_cut

def plot_fluxes(qsr_parm, qsr_parm_err, qsr_image_num, fit_time_var, fit_time_var_err, galFit=True):
      #for k in range(0,len(qsr_parm[:,2])):
      #   print k, 
      #   print qsr_parm[k,2],
      #   print np.log10(qsr_parm[k,2])
      
      figure(figsize=(8,10))
      #ax1 = subplot(411)
      #errorbar(qsr_image_num, fit_time_var, yerr=fit_time_var_err, fmt='k.')
      #y1,y2 = ax1.get_ylim()
      #ylim(y2,y1)
      labels= ['M1', 'S1', 'M2', 'S2', 'G ']
      count = 0
      for k in [2,3,4,5,7]:
          if(not galFit and k==7): continue
          flux, flux_err = qsr_parm[:,k], qsr_parm_err[:,k]
          ax2 = subplot(311)
          mag, mag_err, flux, flux_err, qsr_image_num_cut, fit_time_var_cut, fit_time_var_err_cut = mag_err_and_cut(flux, flux_err, qsr_image_num, fit_time_var, fit_time_var_err)
          #print mag.shape, mag_err.shape, qsr_image_num_cut.shape
          mag_err = np.sqrt(mag_err**2 + fit_time_var_err_cut**2)
          errorbar(qsr_image_num_cut, mag, yerr=mag_err, fmt=',')
          ylim(-8.,-13.)
          ylabel('Instrument Magnitudes')
          xlabel('Image Number')
          grid(True)

          ax3 = subplot(312)
          errorbar(qsr_image_num_cut, mag - fit_time_var_cut + global_zero_point, yerr=mag_err, fmt=',')
          ylim(22.5,17.5)
          ylabel('Magnitudes')
          xlabel('Image Number')
          grid(True)

          ax4 = subplot(313)
          ax4.set_yscale('log')
          h,b = np.histogram( mag - fit_time_var_cut + global_zero_point, bins=linspace(17.,22., 26))
          plot(b[1:], h+1.e-1, drawstyle='steps', lw=3, label=labels[count], alpha=0.7)
          xlim(22.5,17.5)
          grid(True)
          xlabel('Magnitudes')
          ylabel('Counts')
          legend(loc=2, fontsize=14)

          med = np.median(mag - fit_time_var_cut + global_zero_point)
          MAD = 1.4826*np.median(np.abs(mag - fit_time_var_cut + global_zero_point-med))
          print '%s %d\t%1.3f +/- %1.3f'%(labels[count], len(mag), med, MAD)
          count += 1
      subplots_adjust(hspace=0.3, top=0.95)

      mags = []
      mag_errs = []
      if(galFit==False):
          figure()
          count=0
          for k in [2,3,4,5]:
              flux, flux_err = qsr_parm[:,k], qsr_parm_err[:,k]
              mag, mag_err, flux, flux_err, qsr_image_num_cut, fit_time_var_cut, fit_time_var_err_cut = mag_err_and_cut(flux, flux_err, qsr_image_num, fit_time_var, fit_time_var_err)
              mag_err = np.sqrt(mag_err**2 + fit_time_var_err_cut**2)
              ax = subplot(4,1,count+1)
              med = np.median(mag - fit_time_var_cut + global_zero_point)
              errorbar(qsr_image_num_cut, mag - fit_time_var_cut + global_zero_point - med, yerr=mag_err, fmt='k.', ecolor='r', label=labels[count])
              ylim(+0.6,-0.6)
              ylabel('Magnitude - Avg. Mag')
              xlabel('Image Number')
              grid(True)
              legend(loc=2)
              count+=1
          subplots_adjust(hspace=0.3, top=0.95)

###############################################################################################

def gal_S_CCD(image_num, full_time_var, qsr_image_num, qsr_time_var, qsr_time_var_err, qsr_parm, qsr_parm_err, gal_mag = 20.515, ZP=29.425):
    instrument_mag = gal_mag - ZP + full_time_var
    S_CCD = AL.magnitude2flux(instrument_mag) 

    def mag_err_and_cut(_flux, _flux_err, _qsr_image_num, _fit_time_var, _fit_time_var_err):
      cut = _flux_err/_flux<0.3
      flux = _flux[cut]
      flux_err = _flux_err[cut]
      qsr_image_num_cut = _qsr_image_num[cut]
      fit_time_var_cut  = _fit_time_var[cut]
      fit_time_var_err_cut  = _fit_time_var_err[cut]
      mag = AL.flux2magnitude(flux)
      mag_err = AL.fluxErr2magErr(flux, flux_err)
      return mag, mag_err, flux, flux_err, qsr_image_num_cut, fit_time_var_cut, fit_time_var_err_cut

    k=7
    print qsr_parm.shape
    print qsr_parm_err.shape
    flux, flux_err = qsr_parm[:,k], qsr_parm_err[:,k]
    mag, mag_err, flux, flux_err, qsr_image_num_cut, fit_time_var_cut, fit_time_var_err_cut = mag_err_and_cut(flux, flux_err, qsr_image_num, qsr_time_var, qsr_time_var_err)

    figure()
    errorbar(qsr_image_num_cut, flux, flux_err, fmt='o', ms=7.5, mfc='none', mec = 'r', color='r', mew=1, label='Measured')
    semilogy(image_num, S_CCD, 'k.', label = 'Expected')
    np.savez('Lens_Gal_S_CCD.npz', image_num = image_num, S_CCD = S_CCD)

    #S_CCD = AL.magnitude2flux(instrument_mag+0.349) 
    #semilogy(image_num, S_CCD, '.', color='gray', label = 'Measured')
    #S_CCD = AL.magnitude2flux(instrument_mag-0.349) 
    #semilogy(image_num, S_CCD, '.', color='gray', label = 'Measured')
    legend(loc=1)
    xlabel('Image Number')
    ylabel('Instrument Flux')
    grid(True, which='both')



######################################################

def get_qsr_fnms(npz_filename):
    f = np.load(npz_filename)
    fnm_list  = f['fnm_list']
    print f.keys()
    #print fnm_list
    print 'Number of Phot Files:', len(fnm_list)
    qsr_fnms = []
    cut_qsr_fnms = []

    chi_sq_list = []
    max_chi_list = []
    qsr_parm_list = []
    qsr_parm_err_list = []
    qsr_image_num = []
    time_var = []
    time_var_err = []

    cut_chi_sq_list = []
    cut_max_chi_list = []
    cut_qsr_parm_list = []
    cut_qsr_parm_err_list = []
    cut_qsr_image_num = []
    cut_time_var = []
    cut_time_var_err = []
    mjd_time = []
    print 'File Loop'
    count=0
    for fnm in fnm_list:
        #if(count>100): continue
        qfnm = fnm.replace('Phot', 'Quad_1')
        num = f['image_num'][count]
        if(os.path.isfile(qfnm)):
            chi_sq, max_chi, qsr_parms, qsr_parm_err = qsr_file_quality_parameters(qfnm)
            chi_sq_list.append(chi_sq)
            max_chi_list.append(max_chi)
            qsr_parm_list.append(qsr_parms)
            qsr_parm_err_list.append(qsr_parm_err)
            qsr_fnms.append(qfnm)
            qsr_image_num.append(num)
            time_var.append(f['fit_time_var'][count])
            time_var_err.append(f['fit_time_var_err'][count])
            #print qsr_fnms[-1]
            if(    chi_sq < chi_sq_cut 
               and max_chi < max_chi_cut
               and qsr_parm_err[2]/qsr_parms[2] < frac_unc_M1_M2
               and qsr_parm_err[3]/qsr_parms[3] < frac_unc_M1_M2
               and qsr_parm_err[5]/qsr_parms[5] < frac_unc_M1_M2
               and qsr_parm_err[4]/qsr_parms[4] < frac_unc_M1_M2):            
                cut_qsr_fnms.append(qfnm)
                cut_chi_sq_list.append(chi_sq)
                cut_max_chi_list.append(max_chi)
                cut_qsr_parm_list.append(qsr_parms)
                cut_qsr_parm_err_list.append(qsr_parm_err)
                cut_qsr_image_num.append(num)
                cut_time_var.append(f['fit_time_var'][count])
                cut_time_var_err.append(f['fit_time_var_err'][count])
        count+=1

    print 'Finished File Loop'
    #def gal_S_CCD(       image_num, f   ull_time_var,             qsr_image_num,         qsr_time_var,           qsr_time_var_err,      qsr_parm, qsr_parm_err, gal_mag = 20.515, ZP=29.425):

    gal_S_CCD(np.array(qsr_image_num), np.array(time_var), np.array(cut_qsr_image_num), np.array(cut_time_var), np.array(cut_time_var_err), np.array(cut_qsr_parm_list), np.array(cut_qsr_parm_err_list), gal_mag = 20.515, ZP=global_zero_point)

    plot_qsr_chi(np.array(chi_sq_list), np.array(max_chi_list))
    plot_parms(np.array(qsr_parm_list), np.array(qsr_parm_err_list))
    plot_fluxes(np.array(cut_qsr_parm_list), np.array(cut_qsr_parm_err_list), np.array(cut_qsr_image_num), np.array(cut_time_var), np.array(cut_time_var_err))
    print 'Number of Phot Files:         ', len(fnm_list)
    print 'Number of Qsr Files:          ',  len(qsr_fnms)
    print 'Number of Surviving Qsr Files:',  len(cut_qsr_fnms)
    #exit()
    f.close()    
    return cut_qsr_fnms


#######################################################

######################################################

def get_qsr_fnms2(npz_filename, qsr_fnm_list):
    f = np.load(npz_filename)
    fnm_list  = f['fnm_list']
    print f.keys()
    #print fnm_list
    print 'Number of Phot Files:', len(qsr_fnm_list)
    qsr_fnms = []
    found_qsr_fnms = []
    cut_qsr_fnms = []

    chi_sq_list = []
    max_chi_list = []
    qsr_parm_list = []
    qsr_parm_err_list = []
    qsr_image_num = []
    time_var = []
    time_var_err = []

    cut_chi_sq_list = []
    cut_max_chi_list = []
    cut_qsr_parm_list = []
    cut_qsr_parm_err_list = []
    cut_qsr_image_num = []
    cut_time_var = []
    cut_time_var_err = []
    mjd_obs = []
    print 'File Loop'
    count=0
    for fnm in fnm_list:
        #if(count>100): continue
        qfnm = fnm.replace('Phot', 'Quad_2')
        num = int(qfnm.split('/')[-1].split('_')[1])
        print 'qfnm', num, qfnm
        #if(fnm.replace('Phot', 'Quad_1') not in qsr_fnm_list): continue # do not include this if it did not pass the first round of cuts
        if(os.path.isfile(qfnm)):
            chi_sq, max_chi, qsr_parms, qsr_parm_err = qsr_file_quality_parameters(qfnm)
            if('0926' in qfnm):
                print qfnm
                print chi_sq, max_chi
                print qsr_parms
                print qsr_parm_err
            chi_sq_list.append(chi_sq)
            max_chi_list.append(max_chi)
            qsr_parm_list.append(qsr_parms)
            qsr_parm_err_list.append(qsr_parm_err)
            qsr_fnms.append(qfnm)
            qsr_image_num.append(num)
            time_var.append(f['fit_time_var'][count])
            time_var_err.append(f['fit_time_var_err'][count])
            #print qsr_fnms[-1]
            if(    chi_sq < chi_sq_cut 
               and max_chi < max_chi_cut
               and qsr_parm_err[2]/qsr_parms[2] < frac_unc_M1_M2
               and qsr_parm_err[3]/qsr_parms[3] < frac_unc_M1_M2
               and qsr_parm_err[5]/qsr_parms[5] < frac_unc_M1_M2
               and qsr_parm_err[4]/qsr_parms[4] < frac_unc_M1_M2
               and qfnm.replace('Quad_2', 'Quad_1') in qsr_fnm_list):            
                cut_qsr_fnms.append(qfnm)
                cut_chi_sq_list.append(chi_sq)
                cut_max_chi_list.append(max_chi)
                cut_qsr_parm_list.append(qsr_parms)
                cut_qsr_parm_err_list.append(qsr_parm_err)
                cut_qsr_image_num.append(num)
                cut_time_var.append(f['fit_time_var'][count])
                cut_time_var_err.append(f['fit_time_var_err'][count])
                g = np.load(qfnm)
                mjd_obs.append(g['mjd_obs'])
                g.close()

        count+=1
    mjd_obs = np.array(mjd_obs)
    cut_qsr_parm_list = np.array(cut_qsr_parm_list)
    cut_qsr_parm_err_list = np.array(cut_qsr_parm_err_list)
    cut_qsr_image_num = np.array(cut_qsr_image_num)
    cut_time_var = np.array(cut_time_var)
    cut_time_var_err = np.array(cut_time_var_err)
    print 'Finished File Loop'
    #def gal_S_CCD(       image_num, f   ull_time_var,             qsr_image_num,         qsr_time_var,           qsr_time_var_err,      qsr_parm, qsr_parm_err, gal_mag = 20.515, ZP=29.425):

    plot_qsr_chi(np.array(chi_sq_list), np.array(max_chi_list))
    plot_parms(np.array(qsr_parm_list), np.array(qsr_parm_err_list), galFit = False)
    plot_fluxes(np.array(cut_qsr_parm_list), np.array(cut_qsr_parm_err_list), np.array(cut_qsr_image_num), np.array(cut_time_var), np.array(cut_time_var_err),galFit = False)

    mag_0, mag_err_0, flux_0, flux_err_0, _qsr_image_num_cut, fit_time_var_cut, fit_time_var_err_cut = mag_err_and_cut(cut_qsr_parm_list[:,2], cut_qsr_parm_err_list[:,2], cut_qsr_image_num, cut_time_var, cut_time_var_err)
    mag_1, mag_err_1, flux_1, flux_err_1, _qsr_image_num_cut, fit_time_var_cut, fit_time_var_err_cut = mag_err_and_cut(cut_qsr_parm_list[:,3], cut_qsr_parm_err_list[:,3], cut_qsr_image_num, cut_time_var, cut_time_var_err)
    mag_2, mag_err_2, flux_2, flux_err_2, _qsr_image_num_cut, fit_time_var_cut, fit_time_var_err_cut = mag_err_and_cut(cut_qsr_parm_list[:,4], cut_qsr_parm_err_list[:,4], cut_qsr_image_num, cut_time_var, cut_time_var_err)
    mag_3, mag_err_3, flux_3, flux_err_3, _qsr_image_num_cut, fit_time_var_cut, fit_time_var_err_cut = mag_err_and_cut(cut_qsr_parm_list[:,5], cut_qsr_parm_err_list[:,5], cut_qsr_image_num, cut_time_var, cut_time_var_err)

    print mjd_obs.shape, mag_0.shape, mag_1.shape ,mag_2.shape, mag_3.shape, mag_err_0.shape, mag_err_1.shape, mag_err_2.shape, mag_err_3.shape

    np.savez('qsr_light_curves_final.npz', mjd_obs = mjd_obs,
                                           mag_0    = mag_0 - fit_time_var_cut + global_zero_point,
                                           mag_1    = mag_1 - fit_time_var_cut + global_zero_point,
                                           mag_2    = mag_2 - fit_time_var_cut + global_zero_point,
                                           mag_3    = mag_3 - fit_time_var_cut + global_zero_point,
                                           mag_err_0    = np.sqrt(mag_err_0**2 + fit_time_var_err_cut**2),
                                           mag_err_1    = np.sqrt(mag_err_1**2 + fit_time_var_err_cut**2),
                                           mag_err_2    = np.sqrt(mag_err_2**2 + fit_time_var_err_cut**2),
                                           mag_err_3    = np.sqrt(mag_err_3**2 + fit_time_var_err_cut**2))

    print 'Number of Phot Files:         ', len(qsr_fnm_list)
    print 'Number of Qsr Files:          ',  len(qsr_fnms)
    print 'Number of Surviving Qsr Files:',  len(cut_qsr_fnms)
    #exit()
    f.close()    
    return qsr_fnms

#######################################################

absolute_cal('phot_cal.npz')
qsr_fnm_list = get_qsr_fnms('phot_cal.npz')
qsr_fnm_list2 = get_qsr_fnms2('phot_cal.npz', qsr_fnm_list)


show()

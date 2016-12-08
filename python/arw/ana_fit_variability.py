
from pylab import *
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

f = np.load('phot_cal.npz')
image_num = f['image_num']
star_table = ascii.read('sample.csv')
fit_mag = f['fit_mag']
star_index = f['star_index']
star_mag_matrix = f['star_mag_matrix']
star_mag_err_matrix = f['star_mag_err_matrix']
star_mask_matrix = f['star_mask_matrix']
fit_time_var = f['fit_time_var']
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

show()

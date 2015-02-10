from AnalysisLibrary import *
import pylab as plt

plt.rcParams['font.size']=14
plt.rcParams['legend.fontsize']=14
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize']=(10,8)
#ra = 04:38:14.9000 
ra = (4.+38./60.+14.9/60./60.)/24*360 
#dec = -12:17:14.40
dec = -12. - (17./60 +14.4/60./60.)
print 'ra',ra,'dec',dec
#fits_file_name = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0123-e90.fits'
#fits_file_name = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0124-e90.fits'
#fits_file_name = '../data/HE0435_LCOGT/lsc1m009-fl03-20141201-0115-e90.fits'
#fits_file_name = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0123-e90.fits'
fits_file_name1 = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0123-e90_wcs_corrected.fits'
fits_file_name2 = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0124-e90_wcs_corrected.fits'
fits_file_name3 = '../data/HE0435_LCOGT/lsc1m009-fl03-20141201-0115-e90_wcs_corrected.fits'

a1 = FITSmanager(fits_file_name1)
a2 = FITSmanager(fits_file_name2)
a3 = FITSmanager(fits_file_name3)

plt.figure()
ax = plt.subplot(111)
ax.set_yscale('log')
#ax.set_xscale('log')
a1.plot_image_values()
a2.plot_image_values()
a3.plot_image_values()
plt.grid(True)
plt.legend(loc=1)

ra = (4.+38./60.+14.9/60./60.)/24*360 
dec = -12. - (17./60 +14.4/60./60.)

ra1 = (4.+38./60.+14.44/60./60.)/24*360
dec1 = -12. - (16./60 +25.4/60./60.)

ra2 = (4.+38./60.+14.60/60./60.)/24*360
dec2 = -12. - (16./60 +37.0/60./60.)

ra3 = (4.+38./60.+12.97/60./60.)/24*360
dec3 = -12. - (17./60 +51.7/60./60.)

#plt.figure()
a1.plot_image()
a1.plot_image(ra_center=ra,dec_center=dec,rad=0.015)


obj  = SourceImage(a1, ra,  dec,  30)
obj1 = SourceImage(a1, ra1, dec1, 30)
obj2 = SourceImage(a1, ra2, dec2, 30)
obj3 = SourceImage(a1, ra3, dec3, 30)

#obj = a1.image_piece(ra, dec, 30)
#obj1 = a1.image_piece(ra1, dec1, 30)
#obj2 = a1.image_piece(ra2, dec2, 30)
#obj3 = a1.image_piece(ra3, dec3, 30)

#plt.figure()
#plt.imshow(test_obj.image, cmap='gray', interpolation='none')

plt.figure()
plt.subplot(222)
plt.imshow(obj.image, cmap='gray', interpolation='none')
plt.title('ra: %1.4f dec: %1.4f'%(ra,dec))
plt.colorbar()
plt.subplot(223)
plt.imshow(obj1.image, cmap='gray', interpolation='none')
plt.title('ra: %1.4f dec: %1.4f'%(ra1,dec1))
plt.colorbar()
plt.subplot(224)
plt.imshow(obj2.image, cmap='gray', interpolation='none')
plt.title('ra: %1.4f dec: %1.4f'%(ra2,dec2))
plt.colorbar()
plt.subplot(221)
plt.imshow(obj3.image, cmap='gray', interpolation='none')
plt.title('ra: %1.4f dec: %1.4f'%(ra3,dec3))
plt.colorbar()


parms1 = obj1.fit_moffat()
parms2 = obj2.fit_moffat()
parms3 = obj3.fit_moffat()

plt.figure(figsize=(12,10))
plt.subplot(223)
obj1.plot_moffat_residual_2D()
plt.subplot(224)
obj2.plot_moffat_residual_2D()
plt.subplot(221)
obj3.plot_moffat_residual_2D()
plt.subplots_adjust(hspace=0.3)
#plt.show()


plt.figure(figsize=(12,10))
plt.subplot(223)
obj1.plot_radial_profiles()
plt.subplot(224)
obj2.plot_radial_profiles()
plt.subplot(221)
obj3.plot_radial_profiles()
plt.subplots_adjust(hspace=0.3)


plt.figure(figsize=(12,10))
plt.subplot(223)
obj1.plot_radial_profile_residuals()
plt.subplot(224)
obj2.plot_radial_profile_residuals()
plt.subplot(221)
obj3.plot_radial_profile_residuals()
plt.subplots_adjust(hspace=0.3)


plt.figure(figsize=(12,10))
plt.subplot(223)
obj1.plot_moffat_chi()
plt.subplot(224)
obj2.plot_moffat_chi()
plt.subplot(221)
obj3.plot_moffat_chi()
plt.subplots_adjust(hspace=0.3)

plt.show()




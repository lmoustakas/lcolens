from AnalysisLibrary import *
import os
import pylab as plt
from datetime import datetime

plt.rcParams['font.size']=14
plt.rcParams['legend.fontsize']=14
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize']=(10,8)
#ra = 04:38:14.9000 
ra = (4.+38./60.+14.9/60./60.)/24*360 
#dec = -12:17:14.40
dec = -12. - (17./60 +14.4/60./60.)
print 'ra',ra,'dec',dec

#data_dir = '/data2/romerowo/lcogt_data/he045-1223_wcs_corrected/corrected_fits/'
data_dir = '/data2/romerowo/lcogt_data/he045-1223_wcs_corrected/'

# GET FILE NAMES IN DATA DIRECTORY
fnms = os.listdir(data_dir)
print 'Number of files:', len(fnms)
fnms_filtered = []
for f in fnms:
  if(len(f.split('-'))>3):
     fnms_filtered.append(f)
fnms = fnms_filtered
print 'Number of files after filtering:', len(fnms)

#for k in range(0,len(fnms)):
#	print fnms[k].split('-')
# PARSE OUT THE TIME ORDERING INFORMATION
fnms_tm = [fnm.split('-')[2]+'-'+fnm.split('-')[3] for fnm in fnms]
# SORT THE FILE NAME LIST IN TIME ORDER
fnms = [x for (y,x) in sorted(zip(fnms_tm, fnms))]
#for k in range(0,len(fnms)):
#	print fnms[k]
#exit()
print 'Number of files after sorting:', len(fnms)

count=0
fout = open('output.dat','w')
for fnm in fnms: 
  #print fnm
  fits_file_name = data_dir + fnm
  FM = FITSmanager(fits_file_name)
  #a2 = FITSmanager(fits_file_name2)
  #a3 = FITSmanager(fits_file_name3)
  plt.figure(1)
  ax = plt.subplot(111)
  ax.set_yscale('log')
  #ax.set_xscale('log')
  #FM.plot_image_values()
  #date_object = datetime.strptime('Jun 1 2005  1:33PM', '%b %d %Y %I:%M%p')

  date = datetime.strptime(FM.hdulist[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
  #print FM.hdulist[0].header['OBJECT'], FM.hdulist[0].header['DATE-OBS'], date
  FM.histogram_image_values(NBINS=20000)
  #print plt.argmax(FM.hist), FM.bin_edges[plt.argmax(FM.hist)]+FM.BW/2.
  if('rp' in FM.hdulist[0].header['FILTER']):
	  plt.plot([date], [FM.bin_edges[plt.argmax(FM.hist)]+FM.BW/2.], 'ro')
  elif('gp' in FM.hdulist[0].header['FILTER']):
	  plt.plot([date], [FM.bin_edges[plt.argmax(FM.hist)]+FM.BW/2.], 'go')
  else:
	  plt.plot([date], [FM.bin_edges[plt.argmax(FM.hist)]+FM.BW/2.], 'ko')

  out_string = '%d\t'%count
  out_string+= '%s\t'%fits_file_name
  out_string+= '%s\t'%FM.hdulist[0].header['DATE-OBS']
  out_string+= '%1.3e\t'%(FM.bin_edges[plt.argmax(FM.hist)]+FM.BW/2.)
  out_string+= '%s\t'%FM.hdulist[0].header['FILTER']
  print out_string
  out_string+='\n'
  fout.write(out_string)

  del FM
  count+=1
  continue

fout.close()
plt.figure(1)
plt.xlabel('Date')
plt.xticks(rotation=45)
plt.ylabel('Peak of Image Counts Distribution')
plt.show()
exit()
plt.figure()
ax = plt.subplot(111)
ax.set_yscale('log')
#ax.set_xscale('log')
a1.plot_image_values()
#a2.plot_image_values()
#a3.plot_image_values()
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




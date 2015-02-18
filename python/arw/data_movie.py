import matplotlib
matplotlib.use('Agg')
from matplotlib import gridspec
from AnalysisLibrary import *
import os
import pylab as plt
from datetime import datetime
import os
from astropy.io import fits


plt.rcParams['font.size']=15
plt.rcParams['legend.fontsize']=15
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize']=(10,8)
#ra = 04:38:14.9000 
ra = (4.+38./60.+14.9/60./60.)/24*360 
#dec = -12:17:14.40
dec = -12. - (17./60 +14.4/60./60.)
print 'ra',ra,'dec',dec

data_dir = '/data2/romerowo/lcogt_data/he045-1223_wcs_corrected/corrected_fits/'
#data_dir = '/data2/romerowo/lcogt_data/he045-1223/'

# GET FILE NAMES IN DATA DIRECTORY
fnms = os.listdir(data_dir)
print 'Number of files:', len(fnms)
fnms_filtered = []
for f in fnms:
  if(len(f.split('-'))>3):
     fnms_filtered.append(f)
fnms = fnms_filtered
print 'Number of files after filtering:', len(fnms)

print 'READING THE FITS FILE OBSERVATION TIMES FOR ORDERING PURPOSES'
mjd_obs_list = []
fnms_list=[]
max_airmass=0.
min_airmass=1.e10
for i in range(len(fnms)):
    if(i%50==0): print '%d of %d files read'%(i,len(fnms))
    d=fnms[i].split('-')
    #print i,fnms[i]
    fnms_list.append(fnms[i])
    fits_file_name = data_dir + fnms[i]
    hdulist = fits.open(fits_file_name)
    mjd_obs_list.append(float(hdulist[0].header['MJD-OBS']))
    if(float(hdulist[0].header['AIRMASS']) > max_airmass): max_airmass = float(hdulist[0].header['AIRMASS'])
    if(float(hdulist[0].header['AIRMASS']) < min_airmass): min_airmass = float(hdulist[0].header['AIRMASS'])
    hdulist.close()

#print min_airmass
#exit()
#fnms_tm = [fnm.split('-')[2]+'-'+fnm.split('-')[3] for fnm in fnms]
print 'SORTING THE FILE NAME LIST IN TIME ORDER'

fnms = [x for (y,x) in sorted(zip(mjd_obs_list, fnms_list))]
#for k in range(0,len(fnms)):
#	print fnms[k]
#exit()
print 'Number of files after sorting:', len(fnms)
count=0
fout = open('output.dat','w')



# ims is a list of lists, each row is a list of artists to draw in the
# current frame; here we are just animating one artist, the image, in
# each frame
count=0
list0=[]
list1=[]
list2=[]
list3=[]
os.system('rm _tmp*.png')
print 'MAKING IMAGE FIGURES FOR ANIMATION'
max_airmass = 0.5*np.ceil(max_airmass/0.5)
for i in range(len(fnms)):
    d=fnms[i].split('-')
    #gs = plt.GridSpec(1, 2, width_ratios=[3, 1]) 
    print '%d of %d \t %s'%(i,len(fnms), fnms[i])
    fnm = fnms[i]
    fits_file_name = data_dir + fnm
    FM = FITSmanager(fits_file_name)
    #print FM.hdulist[0].header['MJD-OBS']
    mjd = float(FM.hdulist[0].header['MJD-OBS'])
    airmass = float(FM.hdulist[0].header['AIRMASS'])
    #print '\tairmass', airmass
    #if(FM.hdulist[0].header['FILTER']=='rp' and mjd>57007.):
    if(FM.hdulist[0].header['FILTER']=='gp' and mjd>57007.):
	    obj  = SourceImage(FM, ra,  dec,  25)
            min_val = np.min(obj.image)
	    max_val = np.max(obj.image)
	    list0.append(mjd)
	    list1.append(max_val)
	    list2.append(min_val)
	    list3.append(airmass)

	    fig = plt.figure(figsize=(10,5))
	    gs0 = gridspec.GridSpec(1, 2)
	    gs00 = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=gs0[1])
	    ax0 = plt.subplot(gs0[0])
	    #im = plt.imshow(obj.image-min_val, cmap=plt.cm.gray_r, interpolation='None', vmax=np.max(obj.image)-min_val)
	    im = plt.imshow(obj.image-min_val, cmap=plt.cm.gray_r, vmax=np.max(obj.image)-min_val)
            plt.suptitle('%s'%fnm)
	    plt.xlabel('pixel')
	    plt.ylabel('pixel')

	    ax1 = plt.subplot(gs00[0])
	    plt.plot(list0,list1,'r.',label='max')
	    plt.plot([mjd],[max_val],'ro',ms=10)
	    plt.plot(list0,list2,'b.',label='min')
	    plt.plot([mjd],[min_val],'bo',ms=10)
	    plt.xlim(57007.,57011.)
	    plt.ylim(0.,3000.)
	    #plt.ylabel('counts')
	    #plt.title('Image Counts')
	    plt.legend(loc=2, title='Image\nCounts', numpoints=1, frameon=False, borderpad=0., handleheight=0., borderaxespad=0.2, handletextpad=-0.2, labelspacing=0.)
	    ax1.set_xticklabels([])
	    ax1.set_yticks([0,1000,2000,3000])
	    #plt.setp(ax1.get_xticklabels(), visible=False)
	    
	    ax2 = plt.subplot(gs00[1])
	    plt.plot(list0,list3,'b.',label='air\nmass')
	    plt.plot([mjd],[airmass],'bo',ms=10)
	    plt.xlim(57007.,57011.)
	    plt.ylim(1.,max_airmass)
	    #plt.ylabel('counts')
	    #plt.title('Image Counts')
	    plt.legend(loc=6, numpoints=1, frameon=False, borderpad=0., handleheight=0., borderaxespad=-0.1, handletextpad=-0.2, labelspacing=0.)
	    plt.xlabel('MJD')
	    ax2.set_yticks(np.arange(1.,max_airmass+0.5,0.5))
	    plt.subplots_adjust(left=0.06, right=0.95, hspace=0.2, wspace=0.15, bottom=0.12)
	    if(i<10):
		fig.savefig('_tmp_00%d.png'%i, dpi=45)
	    if(i>=10 and i<100):
		fig.savefig('_tmp_0%d.png'%i, dpi=45)
	    if(i>=100):
		fig.savefig('_tmp_%d.png'%i, dpi=45)
	    plt.close()
	    del obj
	    del FM
            count+=1

os.system('convert -delay 15 -loop 0 *.png animation.gif')





import matplotlib
matplotlib.use('Agg')
from AnalysisLibrary import *
import os
import pylab as plt
from datetime import datetime
import os

plt.rcParams['font.size']=14
plt.rcParams['legend.fontsize']=14
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


fig = plt.figure()

# ims is a list of lists, each row is a list of artists to draw in the
# current frame; here we are just animating one artist, the image, in
# each frame
count=0
for i in range(len(fnms)):
    d=fnms[i].split('-')
    if(int(d[2])>=20141216 and i!=302 and i!=304 and i!=346 and i!=348):
	    fig = plt.figure()
	    print '%d of %d \t %s'%(i,len(fnms), fnms[i])
	    fnm = fnms[i]
	    fits_file_name = data_dir + fnm
	    FM = FITSmanager(fits_file_name)
	    obj  = SourceImage(FM, ra,  dec,  30)
	    im = plt.imshow(obj.image, cmap=plt.cm.gray_r)
	    #im = plt.imshow(obj.image, cmap=plt.cm.gray_r, interpolation='None')
            plt.text(1,25,'%s'%fnm,color='blue')
	    plt.xlabel('pixel')
	    plt.ylabel('pixel')
	    fig.savefig('_tmp_%d.png'%i, dpi=50)
	    plt.close()
	    del obj
	    del FM
            count+=1

os.system('convert -delay 30 -loop 0 *.png animation.gif')





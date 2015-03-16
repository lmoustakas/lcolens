import os
from astropy.io import ascii

# SET DATA DIRECTORY
dirc = '/data2/romerowo/lcogt_data/he045-1223_wcs_corrected/'

# LOAD TIME ORDERED DATA WITH BACKGROUND COUNTS AND READ NOISE ESTIMATES
filename_table = ascii.read('t_ord_image_stats.dat')

qsr_reject_list = range(400,419)
#for n in range(nr, len(filename_table['filename'])):
for n in range(0, len(filename_table['filename'])):
  mjd_start = 57008.  
  if(filename_table['mjd'][n]< mjd_start):
	continue
  if(n in qsr_reject_list):
	continue
  #'''
  if(n>40):
	continue
  #'''
  num='%d'%n
  if(n<10):
	num='00%s'%n
  if(n>=10 and n<100):
	num='0%s'%n
  if(n>=100 and n<1000):
	num='%s'%n

  print '\nimage_%s %s'%(num,filename_table['filename'][n])
  os.system('./quadPhot.py -i %s -o image_%s'%(filename_table['filename'][n], num))
  

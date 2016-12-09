import os
import numpy as np
from astropy.io import ascii
import time

import argparse

os.nice(20)

field_inputs_file = '/halo_nobackup/lenssim/romerowo/lcolens/python/arw/LCOGT_HE0435_inputs.txt'
if __name__ == "__main__":
	parser=argparse.ArgumentParser(description='quadPhot_loop routine')
	parser.add_argument("-s","--start", default=0,type=int)
	parser.add_argument("-e","--end",   default=-1,type=int)

	args=parser.parse_args()

	# LOAD TIME ORDERED DATA WITH BACKGROUND COUNTS AND READ NOISE ESTIMATES
	table_data_dirs = []
	table_data_dirs.append('/halo_nobackup/lenssim/romerowo/lcogt_data/he045-1223_wcs_corrected/')
	table_data_dirs.append('/halo_nobackup/lenssim/romerowo/lcogt_data/Dec_2015_100_hours/')

	run_tables = []
	run_tables.append(ascii.read('t_ord_image_stats.dat')) # Dec 2014 data
	run_tables.append(ascii.read('time_ordered_file_list.dat')) # Dec 2015 data

	#filename_table = ascii.read('t_ord_image_stats.dat')
	#filename_table = ascii.read('time_ordered_file_list.dat')

	f = np.load('/halo_nobackup/lenssim/romerowo/lcolens/python/arw/Lens_Gal_S_CCD.npz')
	#print f.keys()
	image_num = f['image_num']
	S_CCD = f['S_CCD']
	#exit()

	#for n in range(nr, len(filename_table['filename'])):
	mjd_start = 57008.  

	run_entries = []
	run_entry_nums = []
	count = -1 # initialize so that the first one will be zero.
	table_id = -1
	for table in run_tables:
	    table_id += 1
	    for k in range(0, len(table['filename'])):
	      count += 1
	      num ='{:04d}'.format(count)
          # PSF 
	      '''
	      fout = open('./inputs/joint_run2/PSF/input.%d'%count, 'w')
	      com  = '/halo_nobackup/lenssim/romerowo/lcolens/python/arw/quadPhot.py -i %s '%table['filename'][k]
	      com += ' '
	      com += '-d %s'%table_data_dirs[table_id]
	      com += ' '
	      com += '-fi %s'%field_inputs_file
	      com += ' '
	      com += '-p 3'
	      com += ' '
          #com += '-sl /halo_nobackup/lenssim/romerowo/lcolens/python/arw/proofStars.txt'
          #com += ' '
	      com += '-o /halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/PSF/image_%s'%num
	      com += ' '
	      com += '>  /halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/PSF/run_%s.log'%num
	      fout.write(com)
	      fout.close()
          '''
          
          # Phot
	      '''
          #./quadPhot.py -L 1 -pf /halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/PSF/image_0395_PSF.npz
	      fout = open('./inputs/joint_run2/Phot/input.%d'%count, 'w')
	      com  = '/halo_nobackup/lenssim/romerowo/lcolens/python/arw/quadPhot.py -L 1 '
	      com += ' '
	      com += '-pf %s'%'/halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/PSF/image_%s_PSF.npz'%num
	      com += ' '
	      com += '-fi %s'%field_inputs_file
	      com += ' '
	      com += '-p 0'
	      com += ' '
	      com += '-o /halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/Phot/image_%s'%num
	      com += ' '
	      com += '>  /halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/Phot/run_%s.log'%num
	      fout.write(com)
	      fout.close()
          '''

          # Quad_1
	      '''
	      #./quadPhot.py -L 1 -pf /halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/PSF/image_0395_PSF.npz
	      fout = open('./inputs/joint_run2/Quad_1/input.%d'%count, 'w')
	      com  = '/halo_nobackup/lenssim/romerowo/lcolens/python/arw/quadPhot.py -L 2 '
	      com += ' '
	      com += '-pf %s'%'/halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/PSF/image_%s_PSF.npz'%num
	      com += ' '
	      com += '-fi %s'%field_inputs_file
	      com += ' '
	      com += '-p 3'
	      com += ' '
	      com += '-o /halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/Quad_1/image_%s'%num
	      com += ' '
	      com += '>  /halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/Quad_1/run_%s.log'%num
	      fout.write(com)
	      fout.close()
          '''

          # Quad_2
          #'''
          #./quadPhot.py -L 1 -pf /halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/PSF/image_0395_PSF.npz
	      if int(num) not in image_num:
	        continue
	      if int(num) in image_num:
	        index = np.where(image_num == int(num))[0][0] 
	      #print '%d\t%1.3f'%(image_num[index], S_CCD[index])
	      fout = open('./inputs/joint_run2/Quad_2/input.%d'%count, 'w')
	      com  = '/halo_nobackup/lenssim/romerowo/lcolens/python/arw/quadPhot.py -L 2 '
	      com += ' '
	      com += '-pf %s'%'/halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/PSF/image_%s_PSF.npz'%num
	      com += ' '
	      com += '-fi %s'%field_inputs_file
	      com += ' '
	      com += '-gf %1.3f'%S_CCD[index]
	      com += ' '
	      com += '-p 3'
	      com += ' '
	      com += '-o /halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/Quad_2/image_%s'%num
	      com += ' '
	      com += '>  /halo_nobackup/lenssim/romerowo/lcolens_outputs/20161123/Quad_2/run_%s.log'%num
	      fout.write(com)
	      fout.close()
          #'''


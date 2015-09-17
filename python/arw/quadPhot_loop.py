import os
from astropy.io import ascii
import time

import argparse

if __name__ == "__main__":
	parser=argparse.ArgumentParser(description='quadPhot_loop routine')
	parser.add_argument("-s","--start", default=0,type=int)
	parser.add_argument("-e","--end",   default=-1,type=int)

	args=parser.parse_args()

	# LOAD TIME ORDERED DATA WITH BACKGROUND COUNTS AND READ NOISE ESTIMATES
	filename_table = ascii.read('t_ord_image_stats.dat')

	qsr_reject_list = []
	for k in range(261,266): # No image
		qsr_reject_list.append(k)
	# no image, cosmic ray peak in the FOV
	qsr_reject_list.append(337) # cosmic ray on image
	qsr_reject_list.append(364) # no image
	qsr_reject_list.append(386) # no image and cosmic ray in FOV
	#for k in range(400,419):
	#	qsr_reject_list.append(k)
	print 'qsr_reject_list', qsr_reject_list
	#exit()
	#for n in range(nr, len(filename_table['filename'])):
	mjd_start = 57008.  

	run_entries = []
	run_entry_nums = []
	for n in range(0, len(filename_table['filename'])):
	  num='%d'%n
	  if(n<10):
		num='00%s'%n
	  if(n>=10 and n<100):
		num='0%s'%n
	  if(n>=100 and n<1000):
		num='%s'%n

	  if(filename_table['mjd'][n]< mjd_start):
		continue

	  print '\nimage_%s %s'%(num, filename_table['filename'][n]) 
	  if(n in qsr_reject_list):
		continue
	  '''
	  if(n>40):
		continue
	  '''
	  run_entries.append(n)
	  run_entry_nums.append(num)

	print len(run_entries)
	#exit()
	for k in range(0,len(run_entries)):
          if(k < args.start): continue
          if(k !=-1 and k >= args.end): continue
	  n = run_entries[k]
	  num = run_entry_nums[k]
	  print k, n, num
	  #print '\nimage_%s %s'%(num,filename_table['filename'][n])
	  os.system('./quadPhot.py -i %s -e 1 -o image_%s > run_%s.log'%(filename_table['filename'][n], num, num))
	  time.sleep(10.)
	  os.system('mv *_%s*.npz /disk4/romerowo/lcolens_outputs/20150916/npzfiles'%num)
	  os.system('mv *_%s*.png /disk4/romerowo/lcolens_outputs/20150916/pngfiles'%num)
	  os.system('mv run_%s.log /disk4/romerowo/lcolens_outputs/20150916/logfiles'%num)
	  

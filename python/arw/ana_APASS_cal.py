'''
ARW 05-JAN-2016
This script extracts the derived zero points from lcogt analysis log files and produces a plot vs. time.
'''

from pylab import *
import numpy as np
import glob
rcParams['figure.facecolor']='white'
rcParams['font.size']=16
def readFile(fnm):
	src_file = ''
	filt=''
	airmass=0.
	mjdobs=0.
	APASS_index = []
	APASS_mag = []
	APASS_mag_err = []
	APASS_gr_col_mag = []
	APASS_gr_col_mag_err = []
	Instrument_mag = []
	Instrument_mag_err = []
	for line in file(fnm):
		if('SOURCE_FILENAME' in line): src_file = line.split('=')[1]
		if('FILTER' in line): filt = line.split('=')[1]
		if('AIRMASS' in line): airmass = float(line.split('=')[1])
		if('MJDOBS' in line): mjdobs = float(line.split('=')[1])
		if(line[0]!='#' and '=' not in line):
			APASS_index.append(int(line.split()[0]))
			APASS_mag.append(float(line.split()[1]))
			APASS_mag_err.append(float(line.split()[2]))
			APASS_gr_col_mag.append(float(line.split()[3]))
			APASS_gr_col_mag_err.append(float(line.split()[4]))
			Instrument_mag.append(float(line.split()[5]))
			Instrument_mag_err.append(float(line.split()[6]))
	return src_file, filt, airmass, mjdobs, np.array(APASS_index), np.array(APASS_mag), np.array(APASS_mag_err), np.array(APASS_gr_col_mag), np.array(APASS_gr_col_mag_err), np.array(Instrument_mag), np.array(Instrument_mag_err) 

files = glob.glob('/disk4/romerowo/lcolens_outputs/20151006/textfiles/*.txt')
files = sort(files)
count=0
figure(1)
figure(2)
for fnm in files:
    print fnm
    count+=1
    src_file, filt, airmass, mjdobs, APASS_index, APASS_mag, APASS_mag_err, APASS_gr_col_mag, APASS_gr_col_mag_err, Instrument_mag, Instrument_mag_err  = readFile(fnm)
    mjdobs-=57008.
    col = 'red'
    if(filt.split()[0]=='rp'): col = 'red'
    if(filt.split()[0]=='gp'): col = 'green'
    print src_file.split()[0], '%1.2f'%airmass, '%1.5f'%mjdobs, col
    figure(1)
    try:
        ax=subplot(211)
        errorbar([mjdobs],[np.median(APASS_mag - Instrument_mag)], yerr=[1.4826*np.median(np.abs(APASS_mag - Instrument_mag - np.median(APASS_mag - Instrument_mag)))], fmt='.', color=col, ms=0.01)
        xlabel('mjd - 57008')
        ylabel('$<$m$_{APASS}$ - m$_{I}$$>$')
        subplot(212,sharex=ax)
        plot([mjdobs],[airmass], '.', color=col)
        xlabel('mjd - 57008')
        ylabel('Air mass')
        ylim(2.4,1.)
    except:
        continue
    '''
    for k in range(0,len(APASS_index)):
	    figure(2)
	    if(APASS_index[k] in range(min(APASS_index),max(APASS_index)+1)):
		    if('' in src_file):
			    print src_file, filt.split(), col
			    subplot(221)
			    #errorbar([mjdobs],[APASS_mag[k] - Instrument_mag[k]], yerr=[np.sqrt(APASS_mag_err[k]**2+Instrument_mag_err[k]**2)], fmt= 'ko')
			    plot([mjdobs],[APASS_mag[k] - Instrument_mag[k]], '.', color=col)
			    subplot(224)
			    #errorbar([airmass],[APASS_mag[k] - Instrument_mag[k]], yerr=[np.sqrt(APASS_mag_err[k]**2+Instrument_mag_err[k]**2)], fmt='ko')
			    plot([APASS_mag[k] - Instrument_mag[k]], [airmass], '.', color=col)
			    subplot(223)
			    plot([mjdobs],[airmass], '.', color=col)
    '''
show()



# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 14:03:23 2015

@author: Chris
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy import wcs
from wcsaxes import WCS
import aplpy
import time
import os
import glob
import os.path


#fits_file_name = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0123-e90.fits'
#fits_file_name = '../data/HE0435_LCOGT/lsc1m004-fl04-20141201-0124-e90.fits'
#fits_file_name = '../data/HE0435_LCOGT/lsc1m009-fl03-20141201-0115-e90.fits'

#Empty lists are created and then printed to ensure data does not carry over
#from previous iteration during testing. 
all_int_time = []
print all_int_time
location = []
print location

#topdir = '../'
#exten = '.fits'
 
def step(ext, dirname, names):
    ext = ext.lower()
 
    for name in names:
        if name.lower().endswith('fits'):
            fits_file_name = (os.path.join(dirname, name))
            hdulist=fits.open(fits_file_name)
            start = hdulist[0].header['UTSTART']
            stop = hdulist[0].header['UTSTOP']
            tel_id = hdulist[0].header['TELID']
            if tel_id == '1m0a':
                tel_id = 'Chile' 
            location.append(tel_id)            
            start_time = float(start[0])*36000 + float(start[1])*3600 + float(start[3])*600 + float(start[4])*60 + float(start[6])*10 + float(start[7]) + float(start[9])*.1 + float(start[10])*.01 + float(start[11])*.001
            stop_time = float(stop[0])*36000 + float(stop[1])*3600 + float(stop[3])*600 + float(stop[4])*60 + float(stop[6])*10 + float(stop[7]) + float(stop[9])*.1 + float(stop[10])*.01 + float(stop[11])*.001
            int_time = stop_time - start_time
            delta_time = abs(int_time - 900)
            all_int_time.append(delta_time)
            #For debuggin and testing purposes, print statements are listed below            
            #print start
            #print stop
            #print tel_id
            #print location
            #print start_time
            #print stop_time
            #print int_time
            #print delta_time
            #print all_int_time
            #print ' '
    
os.path.walk('../', step, 'fits')

print all_int_time
print location
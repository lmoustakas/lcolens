
import numpy as np
import matplotlib
import pylab
import astropy.io.ascii
import astropy.io.fits
import astropy.wcs

class PSFglobals:
    def __init__(self, dataDir, inputFile, outputFileTag, plots):
        self.dataDir = dataDir
        self.inputFile = inputFile
        self.outputFileFlag = outputFileTag
        self.plots = plots

def plot_subimage(FITSdata, ra, dec, Npx):
    bigw=astropy.wcs.WCS(FITSdata.header)
    x,y = bigw.wcs_world2pix(ra, dec ,1)    
    yg, xg = np.mgrid[y-Npx/2:y+Npx/2+1,x-Npx/2:x+Npx/2+1]
    d, r = bigw.wcs_pix2world(xg,yg, 1)
    a = FITSdata.data[y-(Npx-1)/2:y+(Npx-1)/2+1,x-(Npx-1)/2:x+(Npx-1)/2+1]
    a = np.rot90(a)
    a = np.flipud(a)
 
    pylab.figure(figsize=(18,6))
    ax = pylab.subplot(121)
    ax_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    pylab.pcolor(r,d,a)
    pylab.colorbar()
    ax.yaxis.set_major_formatter(ax_formatter)
    ax.xaxis.set_major_formatter(ax_formatter)
    pylab.xticks(rotation=45)
    pylab.ylim(np.min(d),np.max(d))
    pylab.xlim(np.min(r),np.max(r))
    ax2 = pylab.subplot(122)
    pylab.pcolor(r,d,np.log10(a))
    pylab.colorbar()
    pylab.ylim(np.min(d),np.max(d))
    pylab.xlim(np.min(r),np.max(r))
    ax_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax2.yaxis.set_major_formatter(ax_formatter)
    ax2.xaxis.set_major_formatter(ax_formatter)
    pylab.xticks(rotation=45)



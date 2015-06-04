
import numpy as np
import matplotlib
import pylab
import astropy.io.ascii
import astropy.io.fits
import astropy.wcs
import matplotlib.widgets

class PSFglobals:
    def __init__(self, dataDir, inputFile, outputFileTag, plots):
        self.dataDir = dataDir
        self.inputFile = inputFile
        self.outputFileFlag = outputFileTag
        self.plots = plots

def plot_subimage(FITSdata, ra, dec, Npx):
    bigw=astropy.wcs.WCS(FITSdata.header)
    x,y = bigw.wcs_world2pix(ra, dec ,1)    
    xg, yg = np.mgrid[y-Npx/2:y+Npx/2+1,x-Npx/2:x+Npx/2+1]
    r, d = bigw.wcs_pix2world(xg,yg, 1)
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

def button_test():
	freqs = np.arange(2, 20, 3)

	fig, ax = pylab.subplots()
	pylab.subplots_adjust(bottom=0.2)
	t = np.arange(0.0, 1.0, 0.001)
	s = np.sin(2*np.pi*freqs[0]*t)
	l, = pylab.plot(t, s, lw=2)


	class Index:
	    ind = 0
	    def next(self, event):
		self.ind += 1
		i = self.ind % len(freqs)
		ydata = np.sin(2*np.pi*freqs[i]*t)
		l.set_ydata(ydata)
		pylab.draw()

	    def prev(self, event):
		self.ind -= 1
		i = self.ind % len(freqs)
		ydata = np.sin(2*np.pi*freqs[i]*t)
		l.set_ydata(ydata)
		pylab.draw()

	callback = Index()
	axprev = pylab.axes([0.7, 0.05, 0.1, 0.075])
	axnext = pylab.axes([0.81, 0.05, 0.1, 0.075])
	bnext = matplotlib.widgets.Button(axnext, 'Next')
	bnext.on_clicked(callback.next)
	bprev = matplotlib.widgets.Button(axprev, 'Previous')
	bprev.on_clicked(callback.prev)


def browse_APASS_subimages(FITSdata, APASS_table, Npx, filename=''):
    bigw=astropy.wcs.WCS(FITSdata.header)
    #fig, ax = pylab.subplots()
    pylab.figure(figsize=(9,6))
    ax = pylab.subplot(111)

    pylab.subplots_adjust(bottom=0.2)
    #t = np.arange(0.0, 1.0, 0.001)
    #s = np.sin(2*np.pi*freqs[0]*t)
    #l, = pylab.plot(t, s, lw=2)
    k=0
    x,y = bigw.wcs_world2pix(APASS_table['radeg'][k], APASS_table['decdeg'][k],1)    
    yg, xg = np.mgrid[y-Npx/2:y+Npx/2+1,x-Npx/2:x+Npx/2+1]
    r, d = bigw.wcs_pix2world(xg,yg, 1)
    a = FITSdata.data[y-(Npx-1)/2:y+(Npx-1)/2+1,x-(Npx-1)/2:x+(Npx-1)/2+1]
    a = np.rot90(a)
    a = np.flipud(a)
    #im = ax.pcolor(r,d,np.log10(a))
    im = ax.imshow(np.log10(a), interpolation='none')
    im.set_cmap('jet')
    pylab.colorbar(im)
    pylab.suptitle(filename)
    ttl = ax.text(.1, 1.05, 'APASS 0', transform = ax.transAxes, va='top')
    k=0
    class Index:
        print 'Index'
        ind = 0
        def next(self, event):
   	    print 'NEXT'
            self.ind += 1
	    if(self.ind>=len(APASS_table)-1): self.ind=0
            k = self.ind
	    x,y = bigw.wcs_world2pix(APASS_table['radeg'][k], APASS_table['decdeg'][k],1)    
	    yg, xg = np.mgrid[y-Npx/2:y+Npx/2+1,x-Npx/2:x+Npx/2+1]
	    r, d = bigw.wcs_pix2world(xg,yg, 1)
	    a = FITSdata.data[y-(Npx-1)/2:y+(Npx-1)/2+1,x-(Npx-1)/2:x+(Npx-1)/2+1]
	    a = np.rot90(a)
	    a = np.flipud(a)
	    im.set_data(np.log10(a))
	    ttl.set_text('APASS %d : %1.3f, %1.3f'%(k,APASS_table['radeg'][k],APASS_table['decdeg'][k]))
            pylab.draw()

        def prev(self, event):
   	    print 'PREV'
            self.ind -= 1
	    if(self.ind<0): self.ind=len(APASS_table)-1
            k = self.ind
	    x,y = bigw.wcs_world2pix(APASS_table['radeg'][k], APASS_table['decdeg'][k],1)    
	    yg, xg = np.mgrid[y-Npx/2:y+Npx/2+1,x-Npx/2:x+Npx/2+1]
	    r, d = bigw.wcs_pix2world(xg,yg, 1)
	    a = FITSdata.data[y-(Npx-1)/2:y+(Npx-1)/2+1,x-(Npx-1)/2:x+(Npx-1)/2+1]
	    a = np.rot90(a)
	    a = np.flipud(a)
	    im.set_data(np.log10(a))
	    ttl.set_text('APASS %d : %1.3f, %1.3f'%(k,APASS_table['radeg'][k],APASS_table['decdeg'][k]))
            pylab.draw()

    callback = Index()
    axprev = pylab.axes([0.7, 0.05, 0.1, 0.075])
    axnext = pylab.axes([0.81, 0.05, 0.1, 0.075])
    bnext = matplotlib.widgets.Button(axnext, 'Next')
    bnext.on_clicked(callback.next)
    bprev = matplotlib.widgets.Button(axprev, 'Previous')
    bprev.on_clicked(callback.prev)  
    pylab.show()
 



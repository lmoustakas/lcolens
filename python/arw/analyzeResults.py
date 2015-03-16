import numpy as np
from pylab import *
import glob

fnames = glob.glob('image_*_results.npz')
fnames = sorted(fnames)
for f in fnames:
	print f

#results = np.load('out_results.npz')
#print npzfile.files
#print npzfile['m1'], npzfile['me1']

m1 = []
me1 = []
for k in range(0,len(fnames)):
	results = np.load(fnames[k])
	m1.append(results['m1'])
	me1.append(results['m1'])

plot(m1, 'ko')
show()

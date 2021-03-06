import numpy as np
from scipy.optimize import root
from scipy.special import gammaincc
from scipy.signal import fftconvolve
import integrate

import pylab as plt


integrate = integrate.integrate()

# Define a Moffat profile at position x,y centered at x0, y0, with parameters alpha and beta
# Position angle is measured in degrees of the semi-major axis East of North (Counterclockwise)
def moffat(x, y, x0, y0, alpha, beta, q = 1.0, posang = 0.0):
    # Normalize to I0 = 1 for now
    # Define the normalization factor
    normfactor = (beta - 1.0) / (np.pi * alpha * alpha)

    # Precompute the sin and cos
    spa = np.sin(posang / 180.0 * np.pi)
    cpa = np.cos(posang / 180.0 * np.pi)

    # Define xprime coordinates in the rotated frame for convenience
    # This uses a standard rotation matrix
    xp = (x - x0) * cpa - (y - y0) * spa
    yp = (x - x0) * spa + (y - y0) * cpa
    # Defined r^2 (no need to take a square root)
    r2 = xp * xp / q / q + yp * yp
    return normfactor * (1.0 + r2 / alpha / alpha) ** (-1.0 * beta)


# Define a pixelized moffat kernel with Simpson's rule integration
def moffat_kernel(x0, y0, alpha, beta, q=1, posang=0.0, nx = 0, ny = 0, subsamp=3):
    # I think this requires nx and ny to be odd
    # That is something to fix in the long run
    fwhm = alpha * 2.0 * (2.0 ** (1.0 / beta) - 1.0) ** 0.5

    # By default make a kernel that is 4 times the fwhm right now and make sure the kernel size is odd
    if nx==0:
        nx = int(np.ceil(4 * fwhm))
        if nx % 2 == 0:
            nx += 1
    if ny==0:
        ny = int(np.ceil(4 * fwhm))
        if ny % 2 == 0:
            ny += 1

    # The slow, but more understandable way
    # Loop over pixels
    # for j, y in enumerate(range(-nx / 2 + 1, nx/2 + 1)):
        # for i, x in enumerate(range(-nx / 2 + 1, nx/2 + 1)):
            # Do the double integral for each pixel
            # The center of the middle pixel is 0, 0
            # Note everything is y, x
    # kern[j, i] = dblint(x0, x - 0.5, x + 0.5, y0, y - 0.5, y + 0.5, alpha, beta, q=q, posang=posang)
    # Subsample the moffat array
    # Say we have 11 pixels in each direction, then 0, 0 is at array index 5, 5
    # If we subsample that by a factor of 10, we want 101 x 101 pixels
    subkern = np.zeros((ny * subsamp - subsamp + 1, nx * subsamp - subsamp + 1 ))
    x = np.linspace(-nx / 2  + 0.5, nx / 2 + 0.5, nx * subsamp + 1  )
    y = np.linspace(-ny / 2  + 0.5, ny / 2 + 0.5, ny * subsamp + 1  )
    x2d, y2d = np.meshgrid(x, y)
    subkern = moffat(x2d, y2d, x0, y0, alpha, beta, q = q, posang = posang)

    # reshape the subkernel so that we can do the integration without a for loop
    kern4d = integrate.make4d(subkern, nx, ny, subsamp)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    return integrate.simp4d(kern4d, dx, dy)

# Define a mixture of Gaussians PSF
def gaussmix(x, y, x0s, y0s, sigmas, norms):
    #Figure out how many Gaussian components you have
    ngauss = len(norms)
    result = None
    for i in range(ngauss):
        this_gauss = norms[i] / (2.0 * np.pi) ** 0.5 / sigmas[i]
        this_gauss *= np.exp(-1.0/2.0/sigmas[i]/sigmas[i]*((x - x0s[i]) ** 2 + (y - y0s[i]) **2))
        if result is not None:
            result += this_gauss
        else:
            result = this_gauss
    return result
    
def gaussmix_kernel(x0s, y0s, sigmas, norms, nx=0, ny=0, subsamp=3):
    maxfwhm = 2.3458 * sigmas.max()
    # By default make a kernel that is 4 times the fwhm right now and make sure the kernel size is odd
    if nx==0:
        nx = int(np.ceil(4 * fwhm))
        if nx % 2 == 0:
            nx += 1
    if ny==0:
        ny = int(np.ceil(4 * fwhm))
        if ny % 2 == 0:
            ny += 1

    subkern = np.zeros((ny * subsamp - subsamp + 1, nx * subsamp - subsamp + 1 ))
    x = np.linspace(-nx / 2  + 0.5, nx / 2 + 0.5, nx * subsamp + 1  )
    y = np.linspace(-ny / 2  + 0.5, ny / 2 + 0.5, ny * subsamp + 1  )
    x2d, y2d = np.meshgrid(x, y)
    subkern = gaussmix(x2d, y2d, x0s, y0s, sigmas, norms)

    # reshape the subkernel so that we can do the integration without a for loop
    kern4d = integrate.make4d(subkern, nx, ny, subsamp)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    return integrate.simp4d(kern4d, dx, dy)

################################################################
# SHAPELETS ####################################################
################################################################

from scipy.special import genlaguerre
from scipy.misc import factorial

def psi_r(n, m, beta, r):
    # m must be positive
    p = (-1.0) ** ((n - m) / 2.0) * beta ** (-1 - m)
    # Make sure we do integer division here so that the factorial is defined
    p *= (factorial((n - m) // 2) / np.pi)**0.5
    p *= (factorial((n + m) // 2)) ** -0.5
    p *= r ** m
    p *= genlaguerre((n - m) // 2, m)(r ** 2.0 / beta ** 2.0)
    p *= np.exp(-r * r / 2 / beta / beta) 
    return p

def num_params(nmax, mmax):
    counter = 0
    for n in range(nmax + 1):
        for m in range(n, -1, -2)[::-1]: # decending order guarantees parity. The [::-1] term puts it back in ascending order
            if m==0:
                counter += 1
            if m>0 and m<=mmax:
                counter += 2 # real and imaginary components
    return counter

def get_nm(nmax, mmax):
    counter = 0
    nm_list = []
    for n in range(nmax + 1):
        for m in range(n, -1, -2)[::-1]: # decending order guarantees parity. The [::-1] term puts it back in ascending order
            if m==0:
                counter += 1
                nm_list.append([n,m,0])
            if m>0 and m<=mmax:
                counter += 2 # real and imaginary components
                nm_list.append([n,m,0])
                nm_list.append([n,m,1])
    return nm_list

def shapelet_image(X, Y, nmax, mmax, x0, y0, beta,  params):
    r = ((X-x0)**2.0 + (Y-y0)**2.0)**0.5
    theta = np.arctan2(Y-y0, X-x0)
    #print 'r.shape, theta.shape', r.shape, theta.shape
    model = np.zeros((len(Y), len(X)))
    counter = 0
    for n in range(nmax + 1):
        for m in range(n, -1, -2)[::-1]: # decending order guarantees parity. The [::-1] term puts it back in ascending order
            #print 'n,m',n,m, params[counter], counter
            if m==0:
                model += params[counter] * psi_r(n, m, beta, r)
                counter += 1
            if m>0 and m<=mmax:
                model += params[counter] *  2.0 * psi_r(n, m, beta, r) * np.cos(m * theta)
                counter += 1
                model += params[counter] *  2.0 * psi_r(n, m, beta, r) * np.sin(m * theta)
                counter +=1
    return model

def shapelet_model(x,y, parameters, nmax, mmax):
    # parse parameters
    x0 = parameters[0]
    y0 = parameters[1]
    beta = parameters[2]
    bkg = parameters[3]
    params = parameters[4:]
    #print x0, y0
    return shapelet_image(x, y, nmax, mmax, x0, y0, beta,  params)

def shapelet_kernel(parameters, nmax, mmax, nx = 0, ny = 0, subsamp=3):
    # I think this requires nx and ny to be odd
    # That is something to fix in the long run

    # By default make a kernel that is 4 times the fwhm right now and make sure the kernel size is odd
    '''
    if nx==0:
        nx = int(np.ceil(4 * beta))
        if nx % 2 == 0:
            nx += 1
    if ny==0:
        ny = int(np.ceil(4 * beta))
        if ny % 2 == 0:
            ny += 1
    '''
    # Subsample the shapelet array
    # Say we have 11 pixels in each direction, then 0, 0 is at array index 5, 5
    # If we subsample that by a factor of 10, we want 101 x 101 pixels
    #subkern = np.zeros((ny * subsamp - subsamp + 1, nx * subsamp - subsamp + 1 ))
    x = np.linspace(-nx / 2  + 0.5, nx / 2 + 0.5, nx * subsamp + 1  )
    y = np.linspace(-ny / 2  + 0.5, ny / 2 + 0.5, ny * subsamp + 1  )
    x2d, y2d = np.meshgrid(x, y)
    subkern = shapelet_model(x2d, y2d, parameters, nmax, mmax)

    #plt.figure()
    #plt.imshow(subkern)
    #print '**', len(x), len(y), min(x), max(x)

    # reshape the subkernel so that we can do the integration without a for loop
    kern4d = integrate.make4d(subkern, nx, ny, subsamp)
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    img_subpix = integrate.simp4d(kern4d, dx, dy)
    return img_subpix



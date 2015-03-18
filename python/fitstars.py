import numpy as np
from psf import moffat_kernel
from scipy.optimize import minimize, curve_fit

def make_substamps(imagedata, x0s, y0s, subnx, subny):
    nstars = len(x0s)
    # Make 3D array of the cutouts of the image
    substamps = np.zeros((nstars, subny, subnx))
    for i in range(nstars):
        xcenter = np.round(x0s[i])
        ycenter = np.round(y0s[i])
        substamps[i] = imagedata[ycenter - subny/2 - 1: ycenter + subny/2, xcenter - subny/2 - 1: xcenter + subny/2]
    return substamps


def fitstars(imagedata, x0s, y0s, readnoise, subnx=31, subny=31):
    # Estimate the sky
    sky = np.median(imagedata)
    # Make 3D array of the cutouts of the image
    substamps = make_substamps(imagedata, x0s, y0s, subnx, subny)
    # Make a 3D array of the uncertainties
    sig_substamps = (substamps + readnoise * readnoise) ** 0.5

    # Make an initial guess for the parameters
    # Assume the psf is ~ in the center of the cutout
    p0 = np.zeros(4 + 4 * len(x0s))

    # Initial alpha
    p0[0] = 2.
    # Initial beta
    p0[1] = 2.
    # Initial q
    p0[2] = 1
    # Initial position angle
    p0[3] = 0.0

    # Estimate the flux in the substamp - sky and assume that is close to the scale
    p0[6::4] = substamps.sum(axis=2).sum(axis=1) - sky * substamps.shape[1] * substamps.shape[2]

    # Put in the sky estimate
    p0[7::4] = sky
    
    # Run scipy.optimize to minimize the chi^2
    results = minimize(star_chi2, p0, args=(substamps, sig_substamps), method='Nelder-Mead')
    # return the best fit parameters
    return results.x

def model_stars(params, nx, ny):
    nstars = (len(params) - 4) / 4
    model = np.zeros((nstars, ny, nx))
 
    # For each star
    for i in range(nstars):
        # Calculate the psf model
        model[i] = moffat_kernel(params[4 + i * 4], params[5 + i * 4], params[0], params[1],
                                 q=params[2], posang=params[3], nx=nx, ny=ny, subsamp=10)
        model[i] *= params[6 + i * 4]
        model[i] += params[7 + i * 4]
    if params[3] > 90.0 or params[3] < -90:
        model[:,:,:] = np.inf
    return model

def star_chi2(params, d, sig):
    # d and sig are 3d arrays with all of substamps
    # params are as follows
    # p[0:4] = alpha, beta, q, posang
    # p[4:] = x0, y0, flux, sky for each star
    nx = d.shape[2]
    ny = d.shape[1]
    
    model = model_stars(params, nx, ny)
    # Find the total chi^2
    chi2 = (d - model) ** 2.0 / sig ** 2.0
    return chi2.sum()

nx = 0
ny = 0
def model_star_wrapper(d, *params):
    result = model_stars(np.array(params), nx, ny)
    return result.ravel()
    
def curvefitstars(imagedata, x0s, y0s, readnoise, subnx=31, subny=31):
    # Estimate the sky
    sky = np.median(imagedata)
    # Make 3D array of the cutouts of the image
    substamps = make_substamps(imagedata, x0s, y0s, subnx, subny)
    # Make a 3D array of the uncertainties
    sig_substamps = (substamps + readnoise * readnoise) ** 0.5

    # Make an initial guess for the parameters
    # Assume the psf is ~ in the center of the cutout
    p0 = np.zeros(4 + 4 * len(x0s))

    # Initial alpha
    p0[0] = 2.
    # Initial beta
    p0[1] = 2.
    # Initial q
    p0[2] = 1
    # Initial position angle
    p0[3] = 0.0

    # Estimate the flux in the substamp - sky and assume that is close to the scale
    p0[6::4] = substamps.sum(axis=2).sum(axis=1) - sky * substamps.shape[1] * substamps.shape[2]

    # Put in the sky estimate
    p0[7::4] = sky
    
    # Run scipy.optimize to minimize the chi^2
    #results = minimize(star_chi2, p0, args=(substamps, sig_substamps), method='Nelder-Mead')
    global nx
    global ny
    nx = subnx
    ny = subny
    popt, pcov = curve_fit(model_star_wrapper, substamps.shape, substamps.ravel(), p0=p0, sigma=sig_substamps.ravel())
    # return the best fit parameters
    return popt, pcov

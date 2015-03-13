import numpy as np
from psf import moffat_kernel
from scipy.optimize import minimize

def fitstars(imagedata, x0s, y0s, readnoise, subnx=31, subny=31):
    # Estimate the sky
    sky = np.median(imagedata)
    # Make 3D array of the cutouts of the image
    substamps = np.zeros((len(x0s), subny, subnx))
    for i in range(len(x0s)):
        xcenter = np.round(x0s[i])
        ycenter = np.round(y0s[i])
        substamps[i] = imagedata[ycenter - subny/2: ycenter + subny/2 + 1, xcenter - subny/2: xcenter + subny/2 + 1]
    # Make a 3D array of the uncertainties
    sig_substamps = (substamps + readnoise * readnoise) ** 0.5

    # Make an initial guess for the parameters
    # Assume the psf is ~ in the center of the cutout
    p0 = np.zeros(4 + 4 * len(x0s))

    # Initial alpha
    p0[0] = 3.
    # Initial beta
    p0[1] = 2.
    # Initial q
    p0[2] = 1
    # Initial position angle
    p0[3] = 0.0

    # Estimate the total flux in the substamp - sky and assume that is close to the psf flux
    p0[6::4] = substamps.sum(axis=2).sum(axis=1) - substamps.shape[1] * substamps.shape[0] * sky

    # Put in the sky estimate
    p0[7::4] = sky

    # Run scipy.optimize to minimize the chi^2
    results = minimize(star_chi2, p0, args=(substamps, sig_substamps), method='Nelder-Mead')
    # return the best fit parameters
    return results.x


def star_chi2(params, d, sig):
    # d and sig are 3d arrays with all of substamps
    # params are as follows
    # p[0:4] = alpha, beta, q, posang
    # p[4:] = x0, y0, flux, sky for each star
    nx = d.shape[2]
    ny = d.shape[1]
    model = np.zeros(d.shape)
    # For each star
    for i in range((len(params) - 4) / 4):
        # Calculate the psf model
        model[i] = moffat_kernel(params[4 + i * 4], params[5 + i * 4], params[0], params[1],
                                 q=params[2], posang=params[3], nx=nx, ny=ny, subsamp=10)
        model[i] *= params[6 + i * 4]
        model[i] += params[7 + i * 4]
    # Find the total chi^2
    chi2 = (d - model) ** 2.0 / sig ** 2.0
    return chi2.sum()


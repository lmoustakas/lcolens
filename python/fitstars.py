import numpy as np
from psf import moffat_kernel, gaussmix_kernel
from scipy.optimize import minimize, curve_fit

class star_modeller:
    def __init__(self, nx, ny, nstars, psftype='moffat', subsamp=3, minimize_type='curvefit'):
        self.nx = nx
        self.ny = ny
        self.nstars = nstars
        self.psftype = psftype
        self.subsamp = subsamp
        self.minimize_type = minimize_type
    
    def model_stars(self, params):
        model = np.zeros((self.nstars, self.ny, self.nx))
        # For each star
        for i in range(self.nstars):
            # Calculate the psf model
            if self.psftype == 'moffat':
                model[i] = moffat_kernel(params[4 + i * 4], params[5 + i * 4], params[0], params[1],
                                         q=params[2], posang=params[3], nx=self.nx, ny=self.ny, subsamp=self.subsamp)
                model[i] *= params[6 + i * 4]
                model[i] += params[7 + i * 4]
            elif self.psftype == 'gaussmix':
                ngauss = (len(params) - self.nstars * 4) / 4
                ngausspars = 4 * ngauss
                model[i] = gaussmix_kernel(params[ngausspars + i * 4] + params[0:ngausspars:4], 
                                        params[ngausspars + 1 + i * 4] + params[1:ngausspars:4], 
                                        params[2:ngausspars:4], params[3:ngausspars:4], 
                                        nx=self.nx, ny=self.ny, subsamp=self.subsamp)

                model[i] *= params[ngausspars + 2 + i * 4]
                model[i] += params[ngausspars + 3 + i * 4]

        if self.psftype == 'moffat' and (params[3] > 90.0 or params[3] < -90):
            # Stop theta from wrapping around
            model[:,:,:] = np.inf
        return model
    
    def model_stars_curvefit(self, xdata, *params):
        return self.model_stars(np.array(params)).ravel()
    
    def star_chi2(self, params, d, sig):
        # d and sig are 3d arrays with all of substamps
        # params are as follows
        # p[0:4] = alpha, beta, q, posang
        # p[4:] = x0, y0, flux, sky for each star
        model = self.model_stars(params, d)
        # Find the total chi^2
        chi2 = (d - model) ** 2.0 / sig ** 2.0
        return chi2.sum()


def make_substamps(imagedata, x0s, y0s, subnx, subny):
    nstars = len(x0s)
    # Make 3D array of the cutouts of the image
    substamps = np.zeros((nstars, subny, subnx))
    for i in range(nstars):
        xcenter = np.round(x0s[i])
        ycenter = np.round(y0s[i])
        substamps[i] = imagedata[ycenter - subny/2: ycenter + subny/2 + 1, xcenter - subny/2: xcenter + subny/2 + 1]
    return substamps


def fitstars(imagedata, x0s, y0s, readnoise, subnx=31, subny=31, psftype='moffat', minimize_type='curvefit',
             subsamp=3, ngauss=3):
    # Estimate the sky
    sky = np.median(imagedata)
    # Make 3D array of the cutouts of the image
    substamps = make_substamps(imagedata, x0s, y0s, subnx, subny)
    # Make a 3D array of the uncertainties
    sig_substamps = (substamps + readnoise * readnoise) ** 0.5

    if psftype == 'moffat':
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
    elif psftype == 'gaussmix':
        ngausspars = 4 * ngauss
        p0 = np.zeros(ngausspars + 4 * len(x0s))
        p0[3:ngausspars:4] = 1.0 / ngauss
        # Initial widths of the gaussians
        p0[2:ngausspars:4] = 5.0 * 2.0**-np.arange(ngauss)
        p0[ngausspars + 2::4] = substamps.sum(axis=2).sum(axis=1) - sky * substamps.shape[1] * substamps.shape[2]
        # Put in the sky estimate
        p0[ngausspars + 3::4] = sky
        print(p0)
    modeller = star_modeller(subnx, subny, len(x0s), psftype=psftype, subsamp=subsamp)
    
    if minimize_type == 'minimize':
        # Run scipy.optimize to minimize the chi^2
        results = minimize(modeller.star_chi2, p0, args=(substamps, sig_substamps), method='Nelder-Mead')
        popt = results.x
    elif minimize_type == 'curvefit':
        popt, _pcov = curve_fit(modeller.model_stars_curvefit, substamps.shape,
                                substamps.ravel(), p0=p0, sigma=sig_substamps.ravel())
    
    # return the best fit parameters
    return popt


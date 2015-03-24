import numpy as np

class integrate:
    def __init__(self):
        self.trapweights = np.zeros((1, 1))
        self.simp2dweights = np.zeros((1, 1))
        self.simp4dweights = np.zeros((1, 1, 1, 1))
        
    # Define a 2D trapezoidal rule for integration
    # If we subsample, this should be good enough and much faster than double quad
    def trap2d(self, d, dx, dy):
        # d is the data 2d array
        # y coordinate is first index
        if self.trapweights.shape != d.shape:
            self.trapweights= 4.0 * np.ones(d.shape)
            self.trapweights[:,0] = 2.0
            self.trapweights[0,:] = 2.0
            self.trapweights[-1, :] = 2.0
            self.trapweights[:, -1] = 2.0
            self.trapweights[0, 0] = 1.0
            self.trapweights[-1, -1] = 1.0
            self.trapweights[0, -1] = 1.0
            self.trapweights[-1, 0] = 1.0
        result = self.trapweights * d
        return 0.25 * dx * dy * result.sum()

    # Simpson's rule for 2d integration
    def simp2d(self, d, dx, dy):
        if self.simp2dweights.shape != d.shape:
            self.simp2dweights = np.ones(d.shape)
            self.simp2dweights[1::2, 1:-1:2] = 16
            self.simp2dweights[1::2, 2:-1:2] = 8
            self.simp2dweights[2::2, 1:-1:2] = 8
            self.simp2dweights[2::2, 2:-1:2] = 4
            self.simp2dweights[0, 1::2] = 4
            self.simp2dweights[0, 2::2] = 2
            self.simp2dweights[-1, 1::2] = 4
            self.simp2dweights[-1, 2::2] = 2
            self.simp2dweights[1::2, 0] = 4
            self.simp2dweights[2::2, 0] = 2
            self.simp2dweights[1::2, -1] = 4
            self.simp2dweights[2::2, -1] = 2
            self.simp2dweights[0, 0]  = 1
            self.simp2dweights[0, -1]  = 1
            self.simp2dweights[-1, 0]  = 1
            self.simp2dweights[-1, -1]  = 1
        result = d * self.simp2dweights
        return dx * dy / 9.0 * result.sum()

    # Integrate a 4D array into a 2D image using Simpson's Rule

    def simp4d(self, d, dx, dy):
        if self.simp4dweights.shape != d.shape:
            self.simp4dweights = np.ones(d.shape[2:])
            self.simp4dweights[1::2, 1:-1:2] = 16
            self.simp4dweights[1::2, 2:-1:2] = 8
            self.simp4dweights[2::2, 1:-1:2] = 8
            self.simp4dweights[2::2, 2:-1:2] = 4
            self.simp4dweights[0, 1::2] = 4
            self.simp4dweights[0, 2::2] = 2
            self.simp4dweights[-1, 1::2] = 4
            self.simp4dweights[-1, 2::2] = 2
            self.simp4dweights[1::2, 0] = 4
            self.simp4dweights[2::2, 0] = 2
            self.simp4dweights[1::2, -1] = 4
            self.simp4dweights[2::2, -1] = 2
            self.simp4dweights[0, 0]  = 1
            self.simp4dweights[0, -1]  = 1
            self.simp4dweights[-1, 0]  = 1
            self.simp4dweights[-1, -1]  = 1
        result = d
        result[:,:] *= self.simp4dweights
        return dx * dy / 9.0 * result.sum(axis=3).sum(axis=2)

    def make4d(self, d, nx, ny, subsamp):
        # This converts a 2d array into a 4d array that is easy to integrate
        # Time for some Python reshape magic
        # This code is a bear to test and debug
        # Don't ask how long it took me to figure this out
        # It is orders of magnitude faster than a for loop though
        result4d = np.zeros((ny, nx, subsamp + 1, subsamp + 1))
        result4d[:, :, :-1, :-1] = d[:-1, :-1].reshape((ny, subsamp, nx, subsamp)).swapaxes(1, 2)

        result4d[:, :, :-1, -1] = d[:-1, subsamp::subsamp].reshape(ny, subsamp, nx).swapaxes(1, 2)
        result4d[:, :, -1, :-1] = d[subsamp::subsamp,:-1].reshape(ny, nx, subsamp)
        result4d[:, :, -1, -1] = d[subsamp::subsamp,subsamp::subsamp]

        return result4d
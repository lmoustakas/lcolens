import numpy as np

# Define a 2D trapezoidal rule for integration
# If we subsample, this should be good enough and much faster than double quad
def trap2d(d, dx, dy):
    # d is the data 2d array
    # y coordinate is first index
    weights= 4.0 * np.ones(d.shape)
    weights[:,0] = 2.0
    weights[0,:] = 2.0
    weights[-1, :] = 2.0
    weights[:, -1] = 2.0
    weights[0, 0] = 1.0
    weights[-1, -1] = 1.0
    weights[0, -1] = 1.0
    weights[-1, 0] = 1.0
    result = weights * d
    return 0.25 * dx * dy * result.sum()

# Simpson's rule for 2d integration
def simp2d(d, dx, dy):
    weights = np.ones(d.shape)
    weights[1::2, 1:-1:2] = 16
    weights[1::2, 2:-1:2] = 8
    weights[2::2, 1:-1:2] = 8
    weights[2::2, 2:-1:2] = 4
    weights[0, 1::2] = 4
    weights[0, 2::2] = 2
    weights[-1, 1::2] = 4
    weights[-1, 2::2] = 2
    weights[1::2, 0] = 4
    weights[2::2, 0] = 2
    weights[1::2, -1] = 4
    weights[2::2, -1] = 2
    weights[0, 0]  = 1
    weights[0, -1]  = 1
    weights[-1, 0]  = 1
    weights[-1, -1]  = 1
    result = d * weights
    return dx * dy / 9.0 * result.sum()

# Integrate a 4D array into a 2D image using Simpson's Rule

def simp4d(d, dx, dy):
    weights = np.ones(d.shape[2:])
    weights[1::2, 1:-1:2] = 16
    weights[1::2, 2:-1:2] = 8
    weights[2::2, 1:-1:2] = 8
    weights[2::2, 2:-1:2] = 4
    weights[0, 1::2] = 4
    weights[0, 2::2] = 2
    weights[-1, 1::2] = 4
    weights[-1, 2::2] = 2
    weights[1::2, 0] = 4
    weights[2::2, 0] = 2
    weights[1::2, -1] = 4
    weights[2::2, -1] = 2
    weights[0, 0]  = 1
    weights[0, -1]  = 1
    weights[-1, 0]  = 1
    weights[-1, -1]  = 1
    result = d
    result[:,:] *= weights
    return dx * dy / 9.0 * result.sum(axis=3).sum(axis=2)

def make4d(d, nx, ny, subsamp):
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
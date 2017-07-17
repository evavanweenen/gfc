"""
Functions for dealing with Probability Density Functions
"""

import numpy as np

from scipy.stats import multivariate_normal, norm

from .mapping import pmap

class GFC_Exception(Exception):
    pass

def entropy(x, f):
    w = np.diff(x)
    f_ = np.copy(f[:-1])
    f_ /= f_.sum()
    H = - np.sum(f_ * np.log(f_ / w))
    return H

def multivariate(mean, covariance, amplitude = 1.):
    pdf = multivariate_normal(mean = mean, cov = covariance)
    pdf.amp = amplitude
    return pdf

def univariate(mean, sigma, amplitude = 1.):
    pdf = norm(loc = mean, scale = sigma)
    pdf.amp = amplitude
    return pdf

def eval_one_PDF(PDF, grid):
    return PDF.amp * PDF.pdf(grid)

def eval_total_PDF(PDFs, grid_boundaries, grid_resolution = 1.0):
    """
    grid_boundaries = ((xmin, xmax), (ymin, ymax), (zmin, zmax))
    """
    x, y, z = grid_boundaries
    x1, x2 = x ; y1, y2 = y ; z1, z2 = z
    X, Y, Z = np.mgrid[x1:x2:grid_resolution, y1:y2:grid_resolution, z1:z2:grid_resolution]
    xyz = np.empty(X.shape + (3,))
    xyz[:, :, :, 0] = X; xyz[:, :, :, 1] = Y; xyz[:, :, :, 2] = Z
    PDFs_evaluated = pmap(eval_one_PDF, PDFs, xyz) # parmap is twice as fast as a listcomp
    total = np.sum(PDFs_evaluated, axis = 0)

    return total

def flatten_to_2D(evaluated, which="all"):
    assert evaluated.ndim == 3, "gaia_fc.pdf.flatten_to_2D: expected 3D input array, instead got {0}D".format(evaluated.ndim)
    if which == "xy":
        result = evaluated.sum(2)
    elif which == "xz":
        result = evaluated.sum(1)
    elif which == "yz":
        result = evaluated.sum(0)
    elif which == "all":
        xy = evaluated.sum(2)
        xz = evaluated.sum(1)
        yz = evaluated.sum(0)
        result = (xy, xz, yz)
    else:
        raise GFC_Exception("gaia_fc.pdf.flatten_to_2D: unrecognised input value for keyword `which`: {0}".format(which))
    return result

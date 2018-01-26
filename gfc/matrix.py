import numpy as np
from numpy import sin, cos, array
from astropy import table
from pygaia.astrometry.constants import auKmYearPerSec as A_v
from pygaia.astrometry.vectorastrometry import phaseSpaceToAstrometry as toastro, normalTriad, astrometryToPhaseSpace as tophase, sphericalToCartesian as tocart
from pygaia.astrometry import coordinates as coords
from .tgas import *

ICRS2gal = coords.CoordinateTransformation(coords.Transformations.ICRS2GAL)
gal2ICRS = coords.CoordinateTransformation(coords.Transformations.GAL2ICRS)

def w(t, plx_col = "parallax", mura_col = "pmra", mudec_col = "pmdec", vrad_col = None, w_col = "w"):
    """
    Calculate w_vector (vrad, v_alpha, v_delta) in ICRS coordinates which is the input for extreme deconvolution.
    Input
        parallax (mas)
        pmra (mas/yr) - aka muphistar, includes cos
        pmdec (mas/yr) - aka mutheta
    Returns
        w-vector (vrad, valpha, vdelta) (km/s, km/s, km/s) 
    """
    w = np.empty((len(t[plx_col]), 3))
    if vrad_col == None:
        w[:,0] = 0
    else:
        w[:,0] = t[vrad_col]
    w[:,1] = A_v / t[plx_col] * t[mura_col] #notice that mura includes the cos
    w[:,2] = A_v / t[plx_col] * t[mudec_col]
    t.add_column(table.Column(w, w_col))
    
def A(t, ra_rad_col = "ra_rad", dec_rad_col = "dec_rad", A_col = "A"):
    """
    Calculate matrix A of normal triad necessary in order to calculate rotation matrix R.
    Input
        ra (rad)
        dec (rad)
    Returns
        A - normal triad matrix
    """    
    ca = cos(t[ra_rad_col]) ; sa = sin(t[ra_rad_col]) ; cd = cos(t[dec_rad_col]) ; sd = sin(t[dec_rad_col])
    A = np.empty((len(t[ra_rad_col]), 3, 3))
    A[:,0,0] = ca * cd
    A[:,0,1] = -sa
    A[:,0,2] = -ca * sd
    A[:,1,0] = sa * cd
    A[:,1,1] = ca
    A[:,1,2] = -sa * sd
    A[:,2,0] = sd
    A[:,2,1] = 0.
    A[:,2,2] = cd
    #A[:,:,0] = [ca * cd, sa * cd, sd]
    #A[:,:,1] = [-sa, ca, 0.]
    #A[:,:,2] = [-ca * sd, -sa * sd, cd]
    #A[:,:,1], A[:,:,2], A[:,:,0] = normalTriad(t[ra_rad_col], t[dec_rad_col])
    t.add_column(table.Column(A, A_col))
    
def R(t, A_col = "A", R_col = "R", R_inv_col = "R^-1"):
    """
    Calculate the projection matrix R.
    Input
        A - normal triad matrix
    Returns
        R - rotation matrix
        R^-1 - inverse rotation matrix
    """
    R_inv = ICRS2gal.rotationMatrix.dot(t[A_col]).swapaxes(0,1) 
    t.add_column(table.Column(R_inv, R_inv_col))
    R = np.linalg.inv(R_inv)
    t.add_column(table.Column(R, name = R_col))

def C(t, ra_col = "ra", dec_col = "dec", plx_col = "parallax", mura_col = "pmra", mudec_col = "pmdec", C_col = "C"):
    """
    Compose the convariance matrix C of astrometric observables from correlations and errors from data.
    Input
        Errors of ra (rad), dec (rad), parallax (mas), pmra (mas/yr), pmdec (mas/yr)
        Correlations between all above named observables (unitless)
    Returns
        C - covariance matrix of astrometric observables
    """
    C = np.empty((len(t[ra_col]), 5, 5))
    var = [ra_col+"_rad", dec_col+"_rad", plx_col, mura_col, mudec_col]
    covar = [ra_col, dec_col, plx_col, mura_col, mudec_col]
    for i in range(5):
        for j in range(5):
            if i == j:
                C[:,i,j] = t[var[i]+"_error"] * t[var[j]+"_error"]
            elif i < j:
                C[:,i,j] = t[var[i]+"_error"] * t[var[j]+"_error"] * t[covar[i]+"_"+covar[j]+"_corr"]
            else:
                C[:,i,j] = C[:,j,i]
    t.add_column(table.Column(C, C_col))

def Q(t, plx_col = "parallax", mura_col = "pmra", mudec_col = "pmdec", Q_col = "Q"):
    """
    Calculate Q, matrix necessary to compute covariance matrix S in same ICRS coordinates as w-vector.
    Input
        parallax (mas)
        pmra (mas/yr)
        pmdec (mas/yr)
    Returns
        Q - matrix
    """
    Q = np.zeros((len(t[plx_col]), 3, 5))
    Q[:,1,2] = -A_v/(t[plx_col]**2.) * t[mura_col]
    Q[:,1,3] = A_v/t[plx_col]
    Q[:,2,2] = -A_v/(t[plx_col]**2.) * t[mudec_col]
    Q[:,2,4] = A_v/t[plx_col]
    t.add_column(table.Column(Q, Q_col))

def S(t, C_col = "C", Q_col = "Q", S_col = "S"):
    """
    Calculate covariance matrix S in same ICRS coordinates as w-vector.
    S = Q*C*Q.T
    Input
        Q - matrix in order to calculate new covariance matrix
        C - covariance matrix of astrometric observables
    Returns
        S - covariance matrix of w-vector in ICRS coordinates
    """    
    S = np.empty((len(t[C_col]), 3, 3))
    for i in range(len(t[C_col])):
        S[i,:,:] = t[Q_col][i,:,:].dot(t[C_col][i,:,:]).dot(t[Q_col][i,:,:].T)
    t.add_column(table.Column(S, S_col))

def UVW(t, w_col = "w", R_inv_col = "R^-1", UVW_col = "UVW"):
    """
    Calculate v-vector, from inverted rotation matrix R^-1 and w-vector
    v = R^-1 * w
    Input
        w-vector - (vrad, valpha, vdelta) in ICRS coordinates
        R^-1 - inverse rotation matrix
    Returns
        v - original uvw velocity vector in galactic coordinates
    """
    UVW = np.empty((len(t[w_col]), 3))
    for i in range(len(t[w_col])):
        UVW[i,:] = t[R_inv_col][i,:,:].dot(t[w_col][i,:])
    t.add_column(table.Column(UVW, UVW_col))

def transformation(t, vrad_col = None):
    w(t, vrad_col = vrad_col)
    A(t)
    R(t)
    C(t)
    Q(t)
    S(t)
    UVW(t)

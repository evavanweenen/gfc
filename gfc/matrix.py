import numpy as np
from numpy import sin, cos, array
from astropy import table
from pygaia.astrometry.constants import auKmYearPerSec as A_v
from pygaia.astrometry.vectorastrometry import phaseSpaceToAstrometry as toastro, normalTriad, astrometryToPhaseSpace as tophase, sphericalToCartesian as tocart
from pygaia.astrometry.coordinates import CoordinateTransformation, Transformations

ICRS2gal = CoordinateTransformation(Transformations.ICRS2GAL)
gal2ICRS = CoordinateTransformation(Transformations.GAL2ICRS)

def w(t, plx_col = "parallax", mura_col = "pmra", mudec_col = "pmdec", vrad_col = None, w_col = "w"):
    """
    Calculate w_vector (vrad, v_alpha, v_delta) in ICRS coordinates which is the input for extreme deconvolution.
    Input
        parallax (mas)
        pmra (mas/yr) - aka muphistar, includes cos
        pmdec (mas/yr) - aka mutheta
    Returns
        w-vector (valpha, vdelta, vrad) (km/s, km/s, km/s) 
    """
    w = np.empty((len(t[plx_col]), 3))
    w[:,0] = A_v / t[plx_col] * t[mura_col] #notice that mura includes the cos
    w[:,1] = A_v / t[plx_col] * t[mudec_col]
    if vrad_col == None:
        w[:,2] = 0
    else:
        w[:,2] = t[vrad_col]
    t.add_column(table.Column(w, w_col))
    
def pqr(t, ra_rad_col = "ra_rad", dec_rad_col = "dec_rad", pqr_col = "pqr"):
    """
    Calculate matrix pqr of normal triad necessary in order to calculate rotation matrix R.
    Input
        ra (rad)
        dec (rad)
    Returns
        pqr - normal triad matrix
    """    
    pqr = np.array(normalTriad(t[ra_rad_col], t[dec_rad_col])).T
    t.add_column(table.Column(pqr, pqr_col))
    
def R(t, pqr_col = "pqr", R_col = "R", R_inv_col = "R_inv"):
    """
    Calculate the projection matrix R.
    Input
        pqr - normal triad matrix
    Returns
        R - rotation matrix
        R^-1 - inverse rotation matrix for ICRS to galactic coordinates
    """
    R_inv = ICRS2gal.rotationMatrix.dot(t[pqr_col]).swapaxes(0,1) 
    t.add_column(table.Column(R_inv, R_inv_col))
    R = np.linalg.inv(R_inv)
    t.add_column(table.Column(R, R_col))

def C(t, ra_col = "ra", dec_col = "dec", plx_col = "parallax", mura_col = "pmra", mudec_col = "pmdec", vrad_col = None, C_col = "C"):
    """
    Compose the convariance matrix C of astrometric observables from correlations and errors from data.
    Input
        Errors of ra (rad), dec (rad), parallax (mas), pmra (mas/yr), pmdec (mas/yr)
        Correlations between all above named observables (unitless)
    Returns
        C - covariance matrix of astrometric observables
    """
    length = 5
    var = [ra_col+"_rad", dec_col+"_rad", plx_col, mura_col, mudec_col]
    covar = [ra_col, dec_col, plx_col, mura_col, mudec_col]
    if vrad_col != None:
        length += 1
        var.append(vrad_col)
        covar.append(vrad_col)
    C = np.empty((len(t[ra_col+'_rad']), length, length))
    for i in range(length):
        for j in range(length):
            if i == j:
                C[:,i,j] = t[var[i]+"_error"] * t[var[j]+"_error"]
            elif i < j:
                C[:,i,j] = t[var[i]+"_error"] * t[var[j]+"_error"] * t[covar[i]+"_"+covar[j]+"_corr"]
            else:
                C[:,i,j] = C[:,j,i]
    t.add_column(table.Column(C, C_col))

def Q(t, plx_col = "parallax", mura_col = "pmra", mudec_col = "pmdec", vrad_col = None, Q_col = "Q"):
    """
    Calculate Q, matrix necessary to compute covariance matrix S in same ICRS coordinates as w-vector.
    Input
        parallax (mas)
        pmra (mas/yr)
        pmdec (mas/yr)
    Returns
        Q - matrix
    """
    length = 5
    if vrad_col != None:
        length += 1
    Q = np.zeros((len(t[plx_col]), 3, length))
    Q[:,0,2] = -A_v/(t[plx_col]**2.) * t[mura_col]
    Q[:,1,2] = -A_v/(t[plx_col]**2.) * t[mudec_col]
    Q[:,0,3] = A_v/t[plx_col]
    Q[:,1,4] = A_v/t[plx_col]
    t.add_column(table.Column(Q, Q_col))

def S(t, C_col = "C", Q_col = "Q", S_col = "S"):
    """
    Calculate covariance matrix S in same ICRS coordinates as w-vector.
    S = Q*C*Q.T
    Input
        Q - Jacobian matrix used to calculate new covariance matrix
        C - covariance matrix of astrometric observables
    Returns
        S - covariance matrix of w-vector in ICRS coordinates
    """    
    C = t[C_col] ; Q = t[Q_col]
    S = np.empty((len(C), 3, 3))
    for i in range(len(C)):
        S[i,:,:] = Q[i,:,:].dot(C[i,:,:]).dot(Q[i,:,:].T)
    t.add_column(table.Column(S, S_col))

def UVW(t, w_col = "w", R_inv_col = "R_inv", UVW_col = "UVW"):
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
        UVW[i] = t[R_inv_col][i,:,:].dot(t[w_col][i])
    t.add_column(table.Column(UVW, UVW_col))

def transformation(t, vrad_col = None):
    w(t, vrad_col = vrad_col) ; print("w")
    pqr(t) ; print("pqr")
    R(t) ; print("R")
    C(t, vrad_col = vrad_col) ; print("C")
    Q(t, vrad_col = vrad_col) ; print("Q")
    S(t) ; print("S")
    UVW(t) ; print("UVW")

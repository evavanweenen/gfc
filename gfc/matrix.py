import numpy as np
from numpy import sin, cos, array
from astropy import table
from pygaia.astrometry.constants import auKmYearPerSec as A_v
from pygaia.astrometry.vectorastrometry import phaseSpaceToAstrometry as toastro, normalTriad, astrometryToPhaseSpace as tophase, sphericalToCartesian as tocart
from pygaia.astrometry import coordinates as coords
ICRS_to_galactic = coords.CoordinateTransformation(coords.Transformations.ICRS2GAL)
gal2ICRS = coords.CoordinateTransformation(coords.Transformations.GAL2ICRS)

def A(alpha, delta):
    """
    A matrix for a single star with ICRS coords (alpha, delta)
    """
    ca = cos(alpha) ; sa = sin(alpha) ; cd = cos(delta) ; sd = sin(delta)
    A = array([[ca * cd, -sa, -ca * sd], [sa * cd, ca, -sa * sd], [sd, 0, cd]])
    return A

def A_many(alphas, deltas):
    """
    A matrix for a many stars with ICRS coords (alpha, delta)
    """
    ca = cos(alphas) ; sa = sin(alphas) ; cd = cos(deltas) ; sd = sin(deltas)
    A = np.empty((len(alphas), 3, 3))
    A[:,0,0] = ca * cd
    A[:,0,1] = -sa
    A[:,0,2] = -ca * sd
    A[:,1,0] = sa * cd
    A[:,1,1] = ca
    A[:,1,2] = -sa * sd
    A[:,2,0] = sd
    A[:,2,1] = 0.
    A[:,2,2] = cd
    return A

def R_inv(A):
    """
    R^-1 matrix for a star with A matrix A
    """
    return ICRS_to_galactic.rotationMatrix.dot(A)

def R_inv_many(As):
    """
    R^-1 matrices for many stars with A matrix A
    """
    return ICRS_to_galactic.rotationMatrix.dot(As).swapaxes(0, 1)

def w(parallax, mu_alpha_star, mu_delta, v_r = 0):
    el1 = v_r
    el2 = (A_v / parallax) * mu_alpha_star
    el3 = (A_v / parallax) * mu_delta
    w = array([[el1], [el2], [el3]])
    return w

def w_many(parallax, mu_alpha_star, mu_delta, v_r = 0):
    w = np.empty((len(parallax), 3))
    w[:,0] = v_r
    w[:,1] = A_v / parallax * mu_alpha_star
    w[:,2] = A_v / parallax * mu_delta
    return w

def Q_star(parallax, mu_alpha_star, mu_delta):
    Q = array([[0., 0., 0., 0, 0.],
    [0., 0., -k/(parallax**2.) * mu_alpha_star, k/parallax, 0.],
    [0., 0., -k/(parallax**2.) * mu_delta, 0., k/parallax]])
    return Q

def Q_star_many(parallax, mu_alpha_star, mu_delta):
    Q = np.empty((len(parallax), 3, 5))
    Q[:,0,:] = 0.
    Q[:,:,:2] = 0.
    Q[:,1,2] = -A_v/(parallax**2.) * mu_alpha_star
    Q[:,1,3] = A_v/parallax
    Q[:,1,4] = 0.
    Q[:,2,2] = -A_v/(parallax**2.) * mu_delta     
    Q[:,2,3] = 0.
    Q[:,2,4] = A_v/parallax
    return Q

def S(C, Q):
    S = Q.dot(C).dot(Q.T)
    return S

def S_many(C, Q):
    S = np.empty((len(C), 3, 3))
    S[:,0,:] = 0.
    S[:,:,0] = 0.
    S[:,1,1] = Q[:,1,2]**2. * C[:,2,2] + 2 * Q[:,1,2] * Q[:,1,3] * C[:,2,3] + Q[:,1,3]**2. * C[:,3,3]
    S[:,2,2] = Q[:,2,2]**2. * C[:,2,2] + 2 * Q[:,2,4] * Q[:,2,2] * C[:,2,4] + Q[:,2,4]**2. * C[:,4,4]
    S[:,1,2] = Q[:,2,2] * Q[:,1,2] * C[:,2,2] + Q[:,2,2] * Q[:,1,3] * C[:,2,3] + Q[:,2,4] * Q[:,1,2] * C[:,2,4] + Q[:,2,4] * Q[:,1,3] * C[:,3,4]
    S[:,2,1] = S[:,1,2]
    return S

def UVW_wR(w, Rinv):
    return Rinv.dot(w)

def UVW_wR_many(w, Rinv):
    UVW = np.empty((len(w), 3))
    UVW[:,0] = w[:,0].T * Rinv[:,0,0] + w[:,1].T * Rinv[:,0,1] + w[:,2].T * Rinv[:,0,2]
    UVW[:,1] = w[:,0].T * Rinv[:,1,0] + w[:,1].T * Rinv[:,1,1] + w[:,2].T * Rinv[:,1,2]
    UVW[:,2] = w[:,0].T * Rinv[:,2,0] + w[:,1].T * Rinv[:,2,1] + w[:,2].T * Rinv[:,2,2]
    return UVW

def add_w(t, plx_col = "parallax", mura_col = "pmra", mudec_col = "pmdec", v_r_col = None, components = True, vector = True):
    if v_r_col is None:
        ws = w_many(t[plx_col], t[mura_col], t[mudec_col])
    else:
        ws = w_many(t[plx_col], t[mura_col], t[mudec_col], v_r = t[v_r_col])
    if components:
        t.add_column(table.Column(data = ws[:, 0], name = "w1", unit = "km / s"))
        t.add_column(table.Column(data = ws[:, 1], name = "w2", unit = "1 / yr"))
        t.add_column(table.Column(data = ws[:, 2], name = "w3", unit = "1 / yr"))
    if vector:
        t.add_column(table.Column(data = ws, name = "w_vec"))

def add_A(t, ra_rad_col = "ra_rad", dec_rad_col = "dec_rad", A_col = "A"):
    As = A_many(t[ra_rad_col], t[dec_rad_col]) 
    t.add_column(table.Column(data = As, name = A_col))

def add_R(t, A_col = "A", R_col = "R", R_inv_col = "R^-1"):
    R_invs = R_inv_many(t[A_col]) 
    t.add_column(table.Column(data = R_invs, name = R_inv_col))
    Rs = np.linalg.inv(R_invs)
    t.add_column(table.Column(data = Rs, name = R_col))

def add_UVW(t, w_vec_col = "w_vec", R_inv_col = "R^-1", components = True, vector = True):
    UVWs = UVW_wR_many(t[w_vec_col], t[R_inv_col])
    if components:
        t.add_column(table.Column(data = UVWs[:, 0], name = "U", unit = "km / s"))
        t.add_column(table.Column(data = UVWs[:, 1], name = "V", unit = "km / s"))
        t.add_column(table.Column(data = UVWs[:, 2], name = "W", unit = "km / s"))
    if vector:
        t.add_column(table.Column(data = UVWs, name = "UVW_vec", unit = "km / s"))

def add_Q(t, plx_col = "parallax", mura_col = "pmra", mudec_col = "pmdec", Q_col = "Q"):
    Qs = Q_star_many(t[plx_col], t[mura_col], t[mudec_col])
    t.add_column(table.Column(data = Qs, name = Q_col))

def add_S(t, C_col = "C", Q_col = "Q", S_col = "S", verbose = True):
    Ss = S_many(t[C_col], t[Q_col])
    t.add_column(table.Column(data = Ss, name = S_col))


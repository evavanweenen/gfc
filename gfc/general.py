"""
General-use functions that did not fit easily into a submodule
These are all imported into the main gfc module, so you do not need to import gfc.general
"""

from __future__ import division

from matplotlib import rcParams

from os import mkdir
from os.path import isdir
from shutil import rmtree

import numpy as np
from numpy import radians, sin, cos, array

import scipy.stats as sp

from scipy.stats import norm, multivariate_normal as m_n, scoreatpercentile

from pygaia.astrometry.vectorastrometry import phaseSpaceToAstrometry as toastro, normalTriad, astrometryToPhaseSpace as tophase, sphericalToCartesian as tocart
from pygaia.astrometry import coordinates as coords

from pygaia.astrometry.constants import auKmYearPerSec as A_v

from extreme_deconvolution import extreme_deconvolution as xd #imported module from bovy

from astropy import table

from time import time

from argparse import ArgumentParser

from astropy import units as u
from astropy.coordinates import Angle

ICRS2gal = coords.CoordinateTransformation(coords.Transformations.ICRS2GAL)
gal2ICRS = coords.CoordinateTransformation(coords.Transformations.GAL2ICRS)

def timestamp():
    """
    Returns a string with the current Unix time.
    """
    return str(int(time()))

def covariant_array(cov):
    """
    Create a 3x3 covariance array from its unique components.
    The input must be of length 6, and will be put into the array as follows:
    cov[0]   cov[1]   cov[2]
    cov[1]   cov[3]   cov[4]
    cov[2]   cov[4]   cov[5]
    """
    return array([[cov[0], cov[1], cov[2]], [cov[1], cov[3], cov[4]], [cov[2], cov[4], cov[5]]])

def XD_arr(t, *cols):
    """
    Convert an astropy table to a numpy array
    
    Parameters
    ----------
    t: astropy.table.table.Table
        Table to take columns from
    *cols: str
        Keys of columns to use
    """
    tarr = t.as_array()
    return np.column_stack([tarr[col] for col in cols])

def XD(y, e, amplitudes, means, covariances, *args, **kwargs):
    """
    Apply extreme deconvolution (wrapper)
    
    Parameters
    ----------
    y: array-like
        Data to model
    e: array-like
        Errors on data
    amplitudes: array-like
        Initial estimates of amplitudes of Gaussian components
    means: array-like
        Initial estimates of amplitudes of Gaussian components
    covariances: array-like
        Initial estimates of covariances of Gaussian components
    *args, **kwargs:
        Additional (keyword) arguments to be passed to extreme_deconvolution.extreme_deconvolution
    
    Returns
    -------
    a, m, c: np.ndarray
        Amplitudes, means, covariances of Gaussian components in the final model
    L: float
        Log-likelihood of the final model
    """
    
    assert y.shape == e.shape[:2], "Dimensions of y and e do not match: {0} and {1}".format(y.shape, e.shape)
    assert len(means) == len(covariances) == len(amplitudes), "lengths of initial amplitudes ({0}), means ({1}) and amplitudes ({2}) do not match.".format(len(amplitudes), len(means), len(covariances))

    a = np.asfortranarray(amplitudes)
    a = a / np.sum(a)
    m = np.asfortranarray(means)
    c = np.asfortranarray(covariances)
    y = np.asfortranarray(y)
    e = np.asfortranarray(e)
    L = xd.extreme_deconvolution(y, e, a, m, c, *args, **kwargs)

    order = a.argsort()[::-1] # descending order of amplitude
    a = a[order] ; m = m[order] ; c = c[order]

    return a, m, c, L

def calc_AIC(lL,comp,N,D,Kmax):
    return -2./N*(N-1.-D-Kmax/2.)*lL+3*comp*D

def calc_MDL(lL,comp,N,D):
    return -lL + comp*D*np.log2(N)/2.

def perform_XD(xarr, xcovar, project, initamps, initmeans, initcovs, wr, Kr, N):
    dim = 3
    logL = np.tile(np.nan, (len(Kr), len(wr)))
    AIC = np.tile(np.nan, (len(Kr), len(wr)))
    MDL = np.tile(np.nan, (len(Kr), len(wr)))
    max_logL = -1e15 ; min_AIC = 1e15 ; min_MDL = 1e15
    for k in Kr:
        print(k)
        for w in range(len(wr)):
            a_xd, m_xd, c_xd, logL[k-1,w] = XD(xarr, xcovar, initamps[k-1], initmeans[k-1], initcovs[k-1], projection=project, w=wr[w]) #extreme deconvolution
            logL[k-1,w] *= N
            AIC[k-1,w] = calc_AIC(logL[k-1,w],k,N,dim,max(Kr))
            MDL[k-1,w] = calc_MDL(logL[k-1,w],k,N,dim)
            if logL[k-1,w] > max_logL:
                max_logL = logL[k-1,w]
                bestK_logL = k ; bestw_logL = wr[w]
                a_logL = a_xd ; m_logL = m_xd ; c_logL = c_xd
            if AIC[k-1,w] < min_AIC:
                min_AIC = AIC[k-1,w]
                bestK_AIC = k ; bestw_AIC = wr[w]
                a_AIC = a_xd ; m_AIC = m_xd ; c_AIC = c_xd
            if MDL[k-1,w] < min_MDL:
                min_MDL = MDL[k-1,w]
                bestK_MDL = k ; bestw_MDL = wr[w]
                a_MDL = a_xd ; m_MDL = m_xd ; c_MDL = c_xd
    a_test = (a_logL, a_AIC, a_MDL)
    m_test = (m_logL, m_AIC, m_MDL)
    c_test = (c_logL, c_AIC, c_MDL)
    bestK = (bestK_logL, bestK_AIC, bestK_MDL)
    bestw = (bestw_logL, bestw_AIC, bestw_MDL)
    if len(a_test[0]) != bestK[0]: print("Error! The amount of components predicted by XD ", len(a_test[0]), " is not equal to ", bestK[0], " the amount of components calculated.")
    return logL, AIC, MDL, a_test, m_test, c_test, bestK, bestw

def add_rad(t, col, unit, errunit):
    print("Transforming to radians for.. ", col)
    col_rad = Angle(t[col], unit).radian 
    col_rad_err = Angle(t[col+"_error"], errunit).radian
    t.add_column(table.Column(col_rad, col + '_rad'))
    t.add_column(table.Column(col_rad_err, col + '_rad_error'))

def r(x, y, z):
    return np.sqrt(x**2 + y**2 + z**2)

def radial_velocity_distribution(PDFs, ra, dec, plx, mura, mudec, covs, vrad_err = 0.):
    vradmin = -750.0
    vradmax = 750.0
    
    p_vrad_marg = []
    p_vrad_cond = []
    p_vrad_cond_temp = []
        
    l, b = ICRS2gal.transformSkyCoordinates(ra, dec)
    pml, pmb = ICRS2gal.transformProperMotions(ra, dec, mura, mudec)
    cov_gal = ICRS2gal.transformCovarianceMatrix(ra, dec, covs)
    
    jacobian = np.array([[-pml*A_v/plx**2., A_v/plx, 0.],
                      [-pmb*A_v/plx**2., 0., A_v/plx]])
    cov_vtan = jacobian.dot(cov_gal[2:5,2:5]).dot(jacobian.T)
    vtan_gal = np.array([pml*A_v/plx, pmb*A_v/plx])
    
    p_gal, q_gal, r_gal = normalTriad(l, b)
    coordtrans = np.vstack((p_gal, q_gal, r_gal))
    
    for N in PDFs:
        #Marginal distribution of vrad
        mean_vrad_marg = r_gal.dot(N.mean)
        sig_vrad_marg = np.sqrt(r_gal.dot(N.cov).dot(r_gal) + vrad_err**2.)
        newPDF_marg = univariate(mean = mean_vrad_marg, sigma = sig_vrad_marg)
        newPDF_marg.amp = N.amp
        p_vrad_marg.append(newPDF_marg)
        # Preparatory calculations for conditional distribution of vrad
        mean_vl = p_gal.dot(N.mean)
        mean_vb = q_gal.dot(N.mean)
        mean_vtan = np.array([mean_vl, mean_vb])
        # Calculate covariance matrix in pqr representation by transforming the xyz representation and add observational errors.
        cov_pqr_gal = coordtrans.dot(N.cov).dot(coordtrans.T)
        cov_pqr_gal[:2,:2] += cov_vtan
        cov_pqr_gal[2,2] += vrad_err**2.
        # Factor the covariance matrix.
        A = cov_pqr_gal[:2,:2]
        B = cov_pqr_gal[2,2]
        C = cov_pqr_gal[2,:2]
        Ainv = np.linalg.inv(A)
        #Conditional distribution of vrad
        mean_vrad_cond = mean_vrad_marg + C.dot(Ainv.dot(vtan_gal - mean_vtan))
        sig_vrad_cond = np.sqrt(B - C.dot(Ainv).dot(C.T))
        newPDF_cond = sp.norm(loc=mean_vrad_cond, scale=sig_vrad_cond)
        newPDF_cond.A = A
        newPDF_cond.mean_vtan = mean_vtan
        prob_vtan_gal = N.amp * sp.multivariate_normal.pdf(vtan_gal, mean = mean_vtan, cov = A)
        newPDF_cond.prob_vtan_gal = prob_vtan_gal
        p_vrad_cond_temp.append(newPDF_cond)
        
    prob_vtan_gal = np.sum(N_cond.prob_vtan_gal for N_cond in p_vrad_cond_temp)
    
    # Conditional distribution of vrad, properly normalized.
    for N, N_cond in zip(PDFs, p_vrad_cond_temp):
        N_cond.q = N.amp * sp.multivariate_normal.pdf(vtan_gal, mean = N_cond.mean_vtan, cov = N_cond.A) / prob_vtan_gal
        p_vrad_cond.append(N_cond)
    
    
    x = np.arange(vradmin, vradmax, 0.0025)     
    p_m = np.array([N_marg.amp * N_marg.pdf(x) for N_marg in p_vrad_marg]).sum(axis=0)
    p_c = np.array([N_cond.q * N_cond.pdf(x) for N_cond in p_vrad_cond]).sum(axis=0)
    

def radial_velocity_distribution_Olivier(PDFs, ra, dec, plx, mura, mudec, C, x = np.arange(-500., 500., 0.0025), v_r_err = 0.):

    p_vr_M = []
    p_vr_C = []
    p_vr_C_temp = []

    l, b = ICRS2gal.transformSkyCoordinates(ra, dec)
    pml, pmb = ICRS_to_galactic.transformProperMotions(ra, dec, mura,mudec)
    C_gal = ICRS_to_galactic.transformCovarianceMatrix(ra, dec, C)
    jacobian = array([[-pml * A_v/plx**2., A_v/plx, 0.   ],
                     [-pmb * A_v/plx**2., 0.   , A_v/plx]])
    cov_v_tan = jacobian.dot(C_gal[2:5,2:5]).dot(jacobian.T)
    v_t_gal = array([pml*A_v/plx, pmb*A_v/plx])
    p_gal, q_gal, r_gal = normalTriad(l, b)
    coordtrans = np.vstack((p_gal, q_gal, r_gal))
    for N in PDFs:
        mean_v_rad_M = r_gal.dot(N.mean)
        sig_v_rad_M  = np.sqrt(r_gal.dot(N.cov).dot(r_gal) + v_r_err**2.)
        newPDF_M = norm(loc = mean_v_rad_M, scale = sig_v_rad_M)
        newPDF_M.amp = N.amp
        p_vr_M.append(newPDF_M)

        mean_v_l = p_gal.dot(N.mean)
        mean_v_b = q_gal.dot(N.mean)
        mean_v_tan = array([mean_v_l, mean_v_b])
        cov_pqr_gal = coordtrans.dot(N.cov).dot(coordtrans.T)
        cov_pqr_gal[:2,:2] += cov_v_tan
        cov_pqr_gal[2,2]   += v_r_err**2.
        A = cov_pqr_gal[:2,:2]
        B = cov_pqr_gal[2, 2]
        C = cov_pqr_gal[2,:2]
        A_inv = np.linalg.inv(A)
        mean_v_rad_C = mean_v_rad_M + C.dot(A_inv.dot(v_t_gal - mean_v_tan))
        sig_v_rad_C = np.sqrt(B - C.dot(A_inv).dot(C.T))
        newPDF_C = norm(loc=mean_v_rad_C, scale=sig_v_rad_C)
        newPDF_C.A = A
        newPDF_C.mean_v_tan = mean_v_tan
        prob_v_tan_gal = N.amp * m_n.pdf(v_t_gal, mean = mean_v_tan, cov = A)
        newPDF_C.prob_v_tan_gal = prob_v_tan_gal
        p_vr_C_temp.append(newPDF_C)

    prob_v_tan_gal = np.sum(N_C.prob_v_tan_gal for N_C in p_vr_C_temp)

    for N, N_C in zip(PDFs, p_vr_C_temp):
        N_C.q = N.amp * m_n.pdf(v_t_gal, mean = N_C.mean_v_tan, cov = N_C.A) / prob_v_tan_gal
        p_vr_C.append(N_C)

    p_m = array([N_M.amp * N_M.pdf(x) for N_M in p_vr_M])
    p_c = array([N_C.q * N_C.pdf(x) for N_C in p_vr_C])
    p_m = p_m.sum(axis=0)
    p_c = p_c.sum(axis=0)

    y = np.vstack((x, p_c, p_m))

    return y

def radial_velocity_best_estimate(ra, dec, plx, mura, mudec, C, PDFs, x = np.arange(-500., 500., 0.0025), v_r_err = 0., which="c"):
    x, conditional, marginalised = radial_velocity_distribution(PDFs, ra, dec, plx, mura, mudec, C, x = x, v_r_err = v_r_err)
    if which.lower() == "m":
        best = x[marginalised.argmax()]
    elif which.lower() == "c":
        best = x[conditional.argmax()]
    elif which.lower() == "b":
        best = (x[marginalised.argmax()], x[conditional.argmax()])
    else:
        raise Exception("gaia_fc.general.radial_velocity_best_estimate: unrecognised value for keyword `which`: `{0}`".format(which))
    return best

def add_column_attribute(t, col, attr, val):
    setattr(t[col], attr, val)

def add_column(t, data, name, **attributes):
    pass_to_Column = {k:v for k,v in attributes.iteritems() if k in ("dtype", "shape", "length", "description", "unit", "format", "meta")}
    set_separately = {k:v for k,v in attributes.iteritems() if k not in pass_to_Column}
    new_column = table.Column(data = data, name = name, **pass_to_Column)
    t.add_column(new_column)
    for key,val in set_separately.iteritems():
        add_column_attribute(t, name, key, val)

def total_velocity(mu_ra, mu_de, parallax, v_r):
    v_ra = k * mu_ra / parallax
    v_de = k * mu_de / parallax
    return r(v_ra, v_de, v_r)

def mean_sigma(x):
    return np.std(x) / np.sqrt(len(x))

def v_vector(U, V, W):
    return array([[U, V, W]]).T
def S2(tab):
    assert len(tab[0]) == 5, "gaia_fc.general.S2: not enough columns in input. Expected 5, got {0}".format(len(tab[0]))
    ras = array([row[0] for row in tab])
    des = array([row[1] for row in tab])
    Us  = array([row[2] for row in tab])
    Vs  = array([row[3] for row in tab])
    Ws  = array([row[4] for row in tab])
    rhats = map_np(rhat, ras, des)
    As = map_np(R_t, rhats)
    vs = map_np(v_vector, Us, Vs, Ws)
    vmean = np.mean(vs, axis=0)
    ps = map_np(np.dot, As, vs)
    pprimes = ps - As.dot(vmean)
    #S_sq = np.mean([r(*p)**2. for p in pprimes])
    S_sq = rse([r(*p)**2. for p in pprimes])
    return S_sq
def sigma_S2(tab):
    assert len(tab[0]) == 5, "gaia_fc.general.sigma_S2: not enough columns in input. Expected 5, got {0}".format(len(tab[0]))
    ras = array([row[0] for row in tab])
    des = array([row[1] for row in tab])
    Us  = array([row[2] for row in tab])
    Vs  = array([row[3] for row in tab])
    Ws  = array([row[4] for row in tab])
    rhats = map_np(rhat, ras, des)
    As = map_np(R_t, rhats)
    vs = map_np(v_vector, Us, Vs, Ws)
    vmean = np.mean(vs, axis=0)
    ps = map_np(np.dot, As, vs)
    pprimes = ps - As.dot(vmean)
    lengths = array([r(*p) for p in pprimes])
    meanp4 = np.mean(lengths**4.)
    meanp22 = np.mean(lengths**2.)**2.
    return np.sqrt((meanp4 - meanp22) / len(tab))

def rse(x):
    """
    Calculate the Robust Scatter Estimate for an array of values (see
    http://dx.doi.org/10.1051/0004-6361/201628714)

    Parameters
    ----------
    x: array-like
        Array of input values.

    Returns
    -------
    rse: float
        RSE value

    The Robust Scatter Estimate (RSE), defined as 0.390152 * (P90-P10),
    where P10 and P90 are the 10th and 90th percentile of the distribution
    of x.
    """
    rse = 0.390152*(scoreatpercentile(x,90)-scoreatpercentile(x,10))
    return rse

def remove_unused_columns(t):
    cols = ["solution_id", "random_index", "astrometric_n_obs_al", "astrometric_n_obs_ac", "astrometric_n_good_obs_al", \
           "astrometric_n_good_obs_ac", "astrometric_n_bad_obs_al", "astrometric_n_bad_obs_ac", "astrometric_delta_q", \
           "astrometric_excess_noise", "astrometric_excess_noise_sig", "astrometric_primary_flag", \
           "astrometric_relegation_factor", "astrometric_weight_al", "astrometric_weight_ac", "astrometric_priors_used", \
           "matched_observations", "duplicated_source", "scan_direction_strength_k1", "scan_direction_strength_k2", \
           "scan_direction_strength_k3", "scan_direction_strength_k4", "scan_direction_mean_k1", "scan_direction_mean_k2", \
           "scan_direction_mean_k3", "scan_direction_mean_k4", "phot_g_n_obs", "phot_g_mean_flux", "phot_g_mean_flux_error"]
    try:
        t.remove_columns(cols)
    except:
        for c in cols:
            try:
                t.remove_column(c)
            except:
                pass

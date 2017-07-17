"""
Extract the covariance matrix for the astrometry from the Hipparcos Catalogue, van Leeuwen 2007 version.
See appendix B of:
    https://ui.adsabs.harvard.edu/#abs/2014A&A...571A..85M/abstract

Anthony Brown September 2016
"""

import pyfits
from os import environ as env
import numpy as np
from gaia_fc.general import radiantomas

def extract_cov(row):
    nu = row["ntr"] - 5
    sigma = np.array([row["e_ra_rad"], row["e_de_rad"], row["e_plx"], row["e_pm_ra"], row["e_pm_de"]])
    nu29 = 2./(9. + nu)
    Q = nu * np.power((np.sqrt(nu29) * row["f2"] + 1. - nu29), 3)
    u = np.sqrt(Q / nu)
    U = np.array([[row["u1"], row["u2"], row["u4"], row["u7"], row["u11"]],
                  [0.,        row["u3"], row["u5"], row["u8"], row["u12"]],
                  [0.,        0.,        row["u6"], row["u9"], row["u13"]],
                  [0.,        0.,        0.,        row["u10"],row["u14"]],
                  [0.,        0.,        0.,        0.,        row["u15"]]])
    for i, s in enumerate(sigma):
        U[i, i] *= u/s

    U_prod = U.T.dot(U)
    U_prod[[0, 1], :] *= radiantomas # RA and Dec rows    -- colon for clarity
    U_prod[:, [0, 1]] *= radiantomas # RA and Dec columns -- double for RARA/RADec/DecDec terms is intentional

    C = np.linalg.inv(U_prod) * u**2. # n.b. second u is the scalar

    return C

def extract_cov_old(args):
    """
    Extract the covariance matrix for the star with the given HIP identifier.

    Parameters
    ----------

    args - command line arguments

    Returns
    -------

    The covariance matrix as a 5x5 matrix.
    """
    radiantomas = 180 * 3600 * 1000 / np.pi
    umat = np.zeros((5,5))
    hdulist=pyfits.open(env.get('HOME')+'/Hipparcos/FVL2007-cat/hip2.fits')
    tbdata=hdulist[1].data
    hip=tbdata.field('HIP')
    index = np.argwhere(hip == args['hip'])[0][0]
    uwIndex = 1
    nobs = tbdata.field('Ntr')[index]
    gof = tbdata.field('F2')[index]
    nu = nobs - 5
    sigma = np.zeros(5)
    sigma[0] = tbdata.field('e_RArad')[index]
    sigma[1] = tbdata.field('e_DErad')[index]
    sigma[2] = tbdata.field('e_Plx')[index]
    sigma[3] = tbdata.field('e_pmRA')[index]
    sigma[4] = tbdata.field('e_pmDE')[index]

    print("Standard errors for HIP {0}:".format(args['hip']))
    for i in range(5):
        "{0:.2f} ".format(sigma[i])
    print()
    Q = nu * np.power((np.sqrt(2 / (9 * nu)) * gof + 1 - (2 / (9 * nu))), 3)
    u = np.sqrt(Q/nu)
    for j in range(5):
        for i in range(j+1):
            umat[i][j] = tbdata.field("UW{0}".format(uwIndex))[index]
            if (i==j):
                umat[i][i] = umat[i][i]*u/sigma[i]
            uwIndex=uwIndex+1

    parvect = np.array([0,0,1,0,0,])

    print()
    print("Ntr and F2")
    print(nobs,gof)
    print()
    temp = np.dot(umat.transpose(),umat)
    for i in range(5):
        for j in range(5):
            if (i==0 or i==1):
                temp[i][j] = temp[i][j]*radiantomas
            if (j==0 or j==1):
                temp[i][j] = temp[i][j]*radiantomas
    covmat = np.linalg.inv(temp)*u*u#np.max([1,u*u])

    print("Covariance matrix:")
    for i in range(5):
        for j in range(5):
            if (i<=1 or j<=1):
                print "{0:+.2e}".format(covmat[i][j])
            else:
                print "{0:+.6f}".format(covmat[i][j])
        print()
    print()

    print("Sqrt of diagonal elements (all in mas)")
    for i in range(5):
        if (i<=1):
            print("{0:+.2f}".format(np.sqrt(covmat[i][i])*radiantomas), end=" ")
        else:
            print("{0:+.2f}".format(np.sqrt(covmat[i][i])), end=" ")
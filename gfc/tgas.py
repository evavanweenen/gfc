"""
Functions specific to TGAS data
"""

import numpy as np
from . import general as gen
from .io import read_votable
from .mapping import *
from astropy import table
import os
from shutil import rmtree

def read_data_as_table(loc = "tgas/tgas_all.vot"):
    t = read_votable(loc, unit_format="vounit")
    ra_rad = np.radians(t["ra"]) ; ra_rad.name = "ra_rad"
    de_rad = np.radians(t["dec"]); de_rad.name = "dec_rad"
    ra_rad_err = np.radians(t["ra_error"]/3.6e6) ; ra_rad_err.name = "ra_rad_error"
    de_rad_err = np.radians(t["dec_error"]/3.6e6) ; de_rad_err.name = "dec_rad_error"
    t.add_columns((ra_rad, ra_rad_err, de_rad, de_rad_err))
    return t

def add_rad(t):
    ra_rad = np.radians(t["ra"]) ; ra_rad.name = "ra_rad"
    de_rad = np.radians(t["dec"]); de_rad.name = "dec_rad"
    ra_rad_err = np.radians(t["ra_error"]/3.6e6) ; ra_rad_err.name = "ra_rad_error"
    de_rad_err = np.radians(t["dec_error"]/3.6e6) ; de_rad_err.name = "dec_rad_error"
    t.add_columns((ra_rad, de_rad, ra_rad_err, de_rad_err))

def add_w(t):
    ws = gen.w_many(t["parallax"], t["pmra"], t["pmdec"])
    t.add_column(table.Column(data = ws[:, 0], name = "w1", unit = "km / s"))
    t.add_column(table.Column(data = ws[:, 1], name = "w2", unit = "1 / yr"))
    t.add_column(table.Column(data = ws[:, 2], name = "w3", unit = "1 / yr"))
    t.add_column(table.Column(data = ws, name = "w_vec"))

def add_A(t):
    As = gen.A_many(t["ra_rad"], t["dec_rad"]) 
    t.add_column(table.Column(data = As, name = "A"))

def add_R(t):
    R_invs = gen.R_inv_many(t["A"]) 
    t.add_column(table.Column(data = R_invs, name = "R^-1"))
    Rs = np.linalg.inv(R_invs)
    t.add_column(table.Column(data = Rs, name = "R"))

def add_UVW(t):
    UVWs = gen.UVW_wR_many(t["w_vec"], t["R^-1"])
    t.add_column(table.Column(data = UVWs[:, 0], name = "U", unit = "km / s"))
    t.add_column(table.Column(data = UVWs[:, 1], name = "V", unit = "km / s"))
    t.add_column(table.Column(data = UVWs[:, 2], name = "W", unit = "km / s"))
    t.add_column(table.Column(data = UVWs, name = "UVW_vec", unit = "km / s"))

def C_many(ra, radec, raplx, rapmra, rapmdec, dec, decplx, decpmra, decpmdec, plx, plxpmra, plxpmdec, pmra, pmrapmdec, pmdec):
    C = np.empty((len(ra), 5, 5))
    C[:,0,0] = ra**2.
    C[:,1,1] = dec**2.
    C[:,2,2] = plx**2.
    C[:,3,3] = pmra**2.
    C[:,4,4] = pmdec**2.
    C[:,0,1] = ra * dec * radec       ; C[:,1,0] = C[:,0,1]
    C[:,0,2] = ra * plx * raplx       ; C[:,2,0] = C[:,0,2]
    C[:,0,3] = ra * pmra * rapmra     ; C[:,3,0] = C[:,0,3]
    C[:,0,4] = ra * pmdec * rapmdec   ; C[:,4,0] = C[:,0,4]
    C[:,1,2] = dec * plx * decplx     ; C[:,2,1] = C[:,1,2]
    C[:,1,3] = dec * pmra * decpmra   ; C[:,3,1] = C[:,1,3]
    C[:,1,4] = dec * pmdec * decpmdec ; C[:,4,1] = C[:,1,4]
    C[:,2,3] = plx * pmra * plxpmra   ; C[:,3,2] = C[:,2,3]
    C[:,2,4] = plx * pmdec * plxpmdec ; C[:,4,2] = C[:,2,4]
    C[:,3,4] = pmra * pmdec * pmrapmdec;C[:,4,3] = C[:,3,4]

    return C

def C_one(ra, radec, raplx, rapmra, rapmdec, dec, decplx, decpmra, decpmdec, plx, plxpmra, plxpmdec, pmra, pmrapmdec, pmdec):
    C = (ra**2, ra*dec*radec, ra*plx*raplx, ra*pmra*rapmra, ra*pmdec*rapmdec, dec**2, dec*plx*decplx, dec*pmra*decpmra, dec*pmdec*decpmdec, plx**2, plx*pmra*plxpmra, plx*pmdec*plxpmdec, pmra**2, pmra*pmdec*pmrapmdec, pmdec**2)
    C = np.array([[C[0], C[1], C[2] , C[3] , C[4]],
                  [C[1], C[5], C[6] , C[7] , C[8]],
                  [C[2], C[6], C[9] , C[10], C[11]],
                  [C[3], C[7], C[10], C[12], C[13]],
                  [C[4], C[8], C[11], C[13], C[14]]])
    return C

def add_C(t):
    Cs = C_many(t["ra_rad_error"], t["ra_dec_corr"], t["ra_parallax_corr"], t["ra_pmra_corr"], t["ra_pmdec_corr"], t["dec_rad_error"], t["dec_parallax_corr"], t["dec_pmra_corr"], t["dec_pmdec_corr"], t["parallax_error"], t["parallax_pmra_corr"], t["parallax_pmdec_corr"], t["pmra_error"], t["pmra_pmdec_corr"], t["pmdec_error"])
    t.add_column(table.Column(data = Cs, name = "C"))

def add_Q(t):
    Qs = gen.Q_star_many(t["parallax"], t["pmra"], t["pmdec"])
    t.add_column(table.Column(data = Qs, name = "Q"))

def add_S(t, split_into = 55000, verbose = True):
    Ss = gen.S_many(t["C"], t["Q"])
    t.add_column(table.Column(data = Ss, name = "S"))

def wrap_v_r(PDFs, row, *args, **kwargs):
    gen.wrap_v_r(PDFs, row["ra_rad"], row["dec_rad"], row["parallax"], row["pmra"], row["pmdec"], row["C"], *args, **kwargs)

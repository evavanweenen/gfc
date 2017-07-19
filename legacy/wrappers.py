"""
Wrappers -- deprecated, to be removed
"""
import numpy as np

from .io import fetch_v_r_from_simbad, flush
from .gplot import RV_distribution
from .general import radial_velocity_distribution

def v_r(PDFs, ra, dec, parallax, mu_ra, mu_dec, C, name = None, saveto = None, x = np.arange(-500., 500., 0.0025), v_r_err_default = 0., *args, **kwargs):
    v_r, v_r_err = None, v_r_err_default
    if name is not None:
        try:
            v_r, v_r_err = fetch_v_r_from_simbad(name)
        except:
            pass

    P = radial_velocity_distribution(PDFs, ra, dec, parallax, mu_ra, mu_dec, C, x = x, v_r_err = v_r_err)
    RV_distribution(P[0], P[1:], saveto = saveto, name = name, real_vr = v_r, labels=["Conditional", "Marginalised"], *args, **kwargs)

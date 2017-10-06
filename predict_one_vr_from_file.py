"""
(Probably) to test how well radial velocities are modelled.
Unfinished.
"""

import gfc
from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("data_file", help = "Table to pull data from")
parser.add_argument("ID", help = "ID of star to predict v_r for")
parser.add_argument("--ID_col", help = "Column in table to identify star by", default = "tycho2_id") 
parser.add_argument("-m", "--mode", help = "Which mode to use (TGAS/Hipparcos/...)", default = "tgas")
parser.add_argument("--min", help = "Minimum v_r to test", type = float, default = -500.0)
parser.add_argument("--max", help = "Maximum v_r to test", type = float, default = 500.0)
parser.add_argument("--step", help = "Step in v_r_ to test", type = float, default = 0.0025)
parser.add_argument("-e", "--error_col", help = "Column with observed v_r error", default = "eHRV")
parser.add_argument("-o", "--observed", help = "Observed value of radial velocity", default = "SIMBAD")
parser.add_argument("--oerr", "--observed_error", help = "Error on observed value of radial velocity", default = 0., type = float)
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

time_before_reading = gfc.time()
t = gfc.io.read_csv(args.data_file)
gfc.add_rad(t)
time_after_reading = gfc.time()
if args.verbose:
    print "Read data in {0:.1f} seconds".format(time_after_reading - time_before_reading)

t = t[t[args.ID_col] == args.ID]

if args.mode.lower() == "tgas":
    C = gfc.tgas.C_many(t["ra_rad_error"], t["ra_dec_corr"], t["ra_parallax_corr"], t["ra_pmra_corr"], t["ra_pmdec_corr"], t["dec_rad_error"], t["dec_parallax_corr"], t["dec_pmra_corr"], t["dec_pmdec_corr"], t["parallax_error"], t["parallax_pmra_corr"], t["parallax_pmdec_corr"], t["pmra_error"], t["pmra_pmdec_corr"], t["pmdec_error"])[0]
else:
    raise NotImplementedError("v_r prediction for mode \"{0}\" is not implemented.".format(args.mode)) 

if args.observed == "SIMBAD":
    v_obs = gfc.io.fetch_v_r_from_simbad(args.ID)

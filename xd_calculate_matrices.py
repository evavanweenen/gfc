"""
Calculate several matrices used in XD.
Run this script to do it once so it doesn't have to be re-done every time (can be quite computationally heavy).
"""

import gfc
from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("data_file", help = "File containing TGAS data")
parser.add_argument("save_folder", help = "Folder in which results will be saved")
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

if args.verbose:
    print "Now reading data..."
time_before_reading = gfc.time()
t = gfc.io.read_csv(args.data_file) # ASSUME CSV -- todo: check
if args.verbose:
    print "{0} lines in data".format(len(t))
    print "Data contains the following keys:"
    print t.keys()

ra_rad = gfc.radians(t["ra"]) ; ra_rad.name = "ra_rad"
de_rad = gfc.radians(t["dec"]); de_rad.name = "dec_rad"
ra_rad_err = gfc.radians(t["ra_error"]/3.6e6) ; ra_rad_err.name = "ra_rad_error"
de_rad_err = gfc.radians(t["dec_error"]/3.6e6) ; de_rad_err.name = "dec_rad_error"
t.add_columns((ra_rad, de_rad, ra_rad_err, de_rad_err))
t.remove_columns(["solution_id", "random_index", "astrometric_n_obs_al", "astrometric_n_obs_ac", "astrometric_n_good_obs_al", \
                  "astrometric_n_good_obs_ac", "astrometric_n_bad_obs_al", "astrometric_n_bad_obs_ac", "astrometric_delta_q", \
                  "astrometric_excess_noise", "astrometric_excess_noise_sig", "astrometric_primary_flag", \
                  "astrometric_relegation_factor", "astrometric_weight_al", "astrometric_weight_ac", "astrometric_priors_used", \
                  "matched_observations", "duplicated_source", "scan_direction_strength_k1", "scan_direction_strength_k2", \
                  "scan_direction_strength_k3", "scan_direction_strength_k4", "scan_direction_mean_k1", "scan_direction_mean_k2", \
                  "scan_direction_mean_k3", "scan_direction_mean_k4", "phot_g_n_obs", "phot_g_mean_flux", "phot_g_mean_flux_error"])
time_after_reading = gfc.time()
if args.verbose:
    print "Read TGAS data in {0:.1f} seconds".format(time_after_reading - time_before_reading)

time_before_matrices = gfc.time()
gfc.tgas.add_w(t)
if args.verbose:
    print "Added w"

gfc.tgas.add_A(t)
if args.verbose:
    print "Added A"

gfc.tgas.add_R(t)
if args.verbose:
    print "Added R and R^-1"

gfc.tgas.add_C(t)
if args.verbose:
    print "Added C"

gfc.tgas.add_Q(t)
if args.verbose:
    print "Added Q"

gfc.tgas.add_S(t)
if args.verbose:
    print "Added S"
time_after_matrices = gfc.time()
if args.verbose:
    print "Calculated matrices in {0:.1f} seconds".format(time_after_matrices - time_before_matrices)

time_before_writing = gfc.time()
gfc.io.write_table_with_separate_arrays(t, saveto_folder = args.save_folder, verbose = args.verbose)
time_after_writing = gfc.time()
if args.verbose:
    print "Written in {0:.1f} seconds".format(time_after_writing - time_before_writing)

import gaia_fc as g

t = g.io.read_csv("tgas/2MASS/plx_small_noMS.csv")
print len(t)
ra_rad = g.radians(t["ra"]) ; ra_rad.name = "ra_rad"
de_rad = g.radians(t["dec"]); de_rad.name = "dec_rad"
ra_rad_err = g.radians(t["ra_error"]/3.6e6) ; ra_rad_err.name = "ra_rad_error"
de_rad_err = g.radians(t["dec_error"]/3.6e6) ; de_rad_err.name = "dec_rad_error"
t.add_columns((ra_rad,de_rad,ra_rad_err,de_rad_err))
t.remove_columns(["solution_id", "random_index", "astrometric_n_obs_al", "astrometric_n_obs_ac", "astrometric_n_good_obs_al", "astrometric_n_good_obs_ac", "astrometric_n_bad_obs_al", "astrometric_n_bad_obs_ac", "astrometric_delta_q", "astrometric_excess_noise", "astrometric_excess_noise_sig", "astrometric_primary_flag", "astrometric_relegation_factor", "astrometric_weight_al", "astrometric_weight_ac", "astrometric_priors_used", "matched_observations", "duplicated_source", "scan_direction_strength_k1", "scan_direction_strength_k2", "scan_direction_strength_k3", "scan_direction_strength_k4", "scan_direction_mean_k1", "scan_direction_mean_k2", "scan_direction_mean_k3", "scan_direction_mean_k4", "phot_g_n_obs", "phot_g_mean_flux", "phot_g_mean_flux_error"])
print "Read TGAS data"

g.tgas.add_w(t)
print "Added w"

g.tgas.add_A(t)
print "Added A"

g.tgas.add_R(t)
print "Added R and R^-1"

g.tgas.add_C(t)
print "Added C"

g.tgas.add_Q(t)
print "Added Q"

g.tgas.add_S(t)
print "Added S"

g.io.write_table_with_separate_arrays(t, saveto_folder="tgas/2MASS/results/")

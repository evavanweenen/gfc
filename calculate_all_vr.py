"""
Predict radial velocities for stars using a model from XD

NOTE: NOT CURRENTLY USABLE
"""

import gaia_fc as g
from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("crossmatch_file", help = "File containing the TGAS/2MASS/RAVE cross-match table")
parser.add_argument("xd_folder", help = "Folder containing results from XD")
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

sqwrange = g.np.arange(0.5, 3.75, 0.25)
wrange = sqwrange**2.

x = len(t) # len(t)
Ps = g.np.tile(g.np.nan, len(t[:x]))
a = g.time()
for w in wrange:
    f = "tgas/2MASS/XD_K_w/w_{w}/K_{K}/".format(w = w, K = K)
    try:
        L_old = g.np.loadtxt(f+"/L_rave")
        #if L_old < -50000:
        #    continue
    except:
        pass
    print K, w
    amps_xd = g.np.load(f+"amplitudes.npy")
    means_xd = g.np.load(f+"means.npy")
    covs_xd = g.np.load(f+"covariances.npy")

    PDFs = map(g.pdf.multivariate, means_xd, covs_xd, amps_xd)

    summary = g.table.Table(names=("ID", "HRV", "eHRV", "best_c", "best_m", "H_c", "H_m", "lower_c", "upper_c", "lower_m", "upper_m", "P"), dtype=[int, float, float, float, float, float, float, float, float, float, float, float])

    for j,row in enumerate(t[:x]):
        RV = g.radial_velocity_distribution(PDFs, row["ra_rad"], row["dec_rad"], row["parallax"], row["pmra"], row["pmdec"], row["C"], x = g.np.arange(-150, 150, 0.1), v_r_err = row["eHRV"])
        #g.np.save("tgas/vr/{0}.npy".format(row["source_id"]), RV[1:])
        #best_c = RV[0][RV[1].argmax()] ; H_c = g.pdf.entropy(RV[0], RV[1])
        #best_m = RV[0][RV[2].argmax()] ; H_m = g.pdf.entropy(RV[0], RV[2])
        #cdf = g.np.cumsum(RV[1:], axis=1) / 10.
        #lower_c = RV[0][g.np.where(cdf[0] >= 0.025)[0][0]]
        #lower_m = RV[0][g.np.where(cdf[1] >= 0.025)[0][0]]
        #upper_c = RV[0][g.np.where(cdf[0] < 0.975)[0][-1]]
        #upper_m = RV[0][g.np.where(cdf[1] < 0.975)[0][-1]]
        Ps[j] = RV[1][g.np.abs(RV[0]-row["HRV"]).argmin()]
        #summary.add_row([row["source_id"], row["HRV"], row["eHRV"], best_c, best_m, H_c, H_m, lower_c, upper_c, lower_m, upper_m, P])
    
    L = g.np.sum(g.np.log(Ps))
    print L
    print >> open(f+"/L_rave", 'w'), L
    print ""
        #summary.write(f+"summary_rave.txt", format="ascii.fixed_width")
print (g.time() - a)/x * len(t) / 3600
#RVs = g.smap_np(g.radial_velocity_distribution, zip(t["ra_rad"], t["dec_rad"], t["parallax"], t["pmra"], t["pmdec"], t["C"]), PDFs) # THIS TAKES VERY VERY LONG
#g.add_column(t, RVs, "v_r")

#t["pml"], t["pmb"] = zip(*[g.ICRS_to_galactic.transformProperMotions(p,th,mp,mth) for p,th,mp,mth in zip(t["ra_rad"], t["dec_rad"], t["pmra"], t["pmdec"])])
#x,y,z,U,V,W = g.tophase(t["l"], t["b"], t["parallax"], t["pml"], t["pmb"], t["r_v"])
#x.name, y.name, z.name, U.name, V.name, W.name = "x", "y", "z", "U", "V", "W"
#t.add_columns((x,y,z,U,V,W))
#print "Calculated individual stars' velocities"
#
##g.io.write_table_without_arrays(t, saveto="tgas/2MASS/TGAS_2MASS_afterXD.dat", exclude_cols=["solution_id", "random_index", "astrometric_n_obs_al", "astrometric_n_obs_ac", "astrometric_n_good_obs_al", "astrometric_n_good_obs_ac", "astrometric_n_bad_obs_al", "astrometric_n_bad_obs_ac", "astrometric_delta_q", "astrometric_excess_noise", "astrometric_excess_noise_sig", "astrometric_primary_flag", "astrometric_relegation_factor", "astrometric_weight_al", "astrometric_weight_ac", "astrometric_priors_used", "matched_observations", "duplicated_source", "scan_direction_strength_k1", "scan_direction_strength_k2", "scan_direction_strength_k3", "scan_direction_strength_k4", "scan_direction_mean_k1", "scan_direction_mean_k2", "scan_direction_mean_k3", "scan_direction_mean_k4", "phot_g_n_obs", "phot_g_mean_flux", "phot_g_mean_flux_error"])
#g.io.write_table_with_separate_arrays(t, saveto_folder="tgas/2MASS/results/", exclude_cols=["solution_id", "random_index", "astrometric_n_obs_al", "astrometric_n_obs_ac", "astrometric_n_good_obs_al", "astrometric_n_good_obs_ac", "astrometric_n_bad_obs_al", "astrometric_n_bad_obs_ac", "astrometric_delta_q", "astrometric_excess_noise", "astrometric_excess_noise_sig", "astrometric_primary_flag", "astrometric_relegation_factor", "astrometric_weight_al", "astrometric_weight_ac", "astrometric_priors_used", "matched_observations", "duplicated_source", "scan_direction_strength_k1", "scan_direction_strength_k2", "scan_direction_strength_k3", "scan_direction_strength_k4", "scan_direction_mean_k1", "scan_direction_mean_k2", "scan_direction_mean_k3", "scan_direction_mean_k4", "phot_g_n_obs", "phot_g_mean_flux", "phot_g_mean_flux_error"])
#print "Written file" 

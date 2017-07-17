"""
Olivier Burggraaff
"""

import gaia_fc as g

t = g.io.read_csv("tgas/tgas_tmass_MS.csv")
ra_rad = g.radians(t["ra"]) ; ra_rad.name = "ra_rad"
de_rad = g.radians(t["dec"]); de_rad.name = "dec_rad"
ra_rad_err = g.radians(t["ra_error"]/3.6e6) ; ra_rad_err.name = "ra_rad_error"
de_rad_err = g.radians(t["dec_error"]/3.6e6) ; de_rad_err.name = "dec_rad_error"
t.add_columns((ra_rad,de_rad,ra_rad_err,de_rad_err))
print "Read TGAS data"
g.gplot.density(t["ra"], t["dec"], saveto="tgas/2MASS/tgas_tmass_sky.png", cb=False, log=False, title="TGAS/2MASS: sky coordinates", xlabel="right ascension", ylabel="declination", cmap="gray")
g.gplot.density(t["l"], t["b"], saveto="tgas/2MASS/tgas_tmass_gal.png", cb=False, log=False, title="TGAS/2MASS: galactic coordinates", xlabel="$\ell$", ylabel="$b$", cmap="gray")
g.gplot.density(t["g_min_ks"], t["g_mag_abs"], saveto="tgas/2MASS/tgas_tmass_HR.png", title="TGAS/2MASS: HR diagram", xlabel="$G - K$ colour", ylabel="Absolute $G$ magnitude", flip="y")
print "Plotted RA/Dec, l/b, HR"

g.tgas.add_w(t)
print "Added w"

g.gplot.density(t["w2"], t["w3"], "tgas/2MASS/w2_w3.png")
print "Plotted w2/w3"

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

warr = g.XD_arr(t, "w1", "w2", "w3")
wcov = g.XD_arr(t, "S") ; wcov[:, 0, 0] = 1e15
proj = g.XD_arr(t, "R")

initial_amps = [24., 23., 23., 9., 9., 7., 4., 2., 1., 0.1] # %

initial_means = [(5., -7., -9.),
                 (-23., -10., -7.),
                 (-13., -30., -8.),
                 (-20., -33., -5.),
                 (9., 4., -6.),
                 (-18., -23., -5.),
                 (-9., -21., -5.),
                 (-40., -19., 0.),
                 (-29., -100., 3.),
                 (2., -100., -10.)]

initial_covs = [(700., -111., -60., 200., 25., 145.),
                (243, 70, 10, 48, 10, 40),
                (1836, 55, -60, 670, -30, 540),
                (350, 165, 110, 230, 110, 135),
                (80, -30, -20, 25, 9, 50),
                (70, -30, -45, 18, 17, 60),
                (9, -7, -6, 13, 8, 14),
                (30, 0.2, 13, 0.6, 0.2, 6),
                (4600, -2400, -330, 4500, -125, 3500),
                (1, 0.5, 0.3, 4, 0.5, 3)]
initial_covs = [g.covariant_array(cov) for cov in initial_covs]

print "Doing XD ...", ; g.io.flush()
amps_xd, means_xd, covs_xd = g.XD(warr, wcov, initial_amps, initial_means, initial_covs, projection = proj, w = 4.)
print "done!"
g.io.save_PDFs(amps_xd, means_xd, covs_xd, "tgas/2MASS/results")

PDFs = map(g.pdf.multivariate, means_xd, covs_xd, amps_xd)
g.gplot.PDFs_multi(PDFs, saveto="tgas/2MASS/PDFs-{0}_{1}.png", title = "Gaussians fitted to TGAS data with XD")
print "Plotted individual PDFs"

evaluated = g.pdf.eval_total_PDF(PDFs, [(-201,201), (-199,199), (-200,200)])
evalxy, evalxz, evalyz = g.pdf.flatten_to_2D(evaluated)
g.gplot.PDFs_gradients(evalxy, vmin=1e-5, extent=(-201,201,-199,199), xlim=(-110,110), ylim=(-120,60), xlabel="$U$ (km/s)", ylabel="$V$ (km/s)", title="TGAS-2MASS cross section $U$ vs $V$ (XD results)", saveto="tgas/2MASS/TGAS_2MASS_XD_PDF_U_V.png")
g.gplot.PDFs_gradients(evalxz, vmin=1e-5, extent=(-201,201,-200,200), xlim=(-110,110), ylim=(-70,70), xlabel="$U$ (km/s)", ylabel="$W$ (km/s)", title="TGAS-2MASS cross section $U$ vs $W$ (XD results)", saveto="tgas/2MASS/TGAS_2MASS_XD_PDF_U_W.png")
g.gplot.PDFs_gradients(evalyz, vmin=1e-5, extent=(-199,199,-200,200), xlim=(-110,110), ylim=(-70,70), xlabel="$V$ (km/s)", ylabel="$W$ (km/s)", title="TGAS-2MASS cross section $V$ vs $W$ (XD results)", saveto="tgas/2MASS/TGAS_2MASS_XD_PDF_V_W.png")
print "Plotted PDF gradients"

RVs = g.smap_np(g.radial_velocity_best_estimate, zip(t["ra_rad"], t["dec_rad"], t["parallax"], t["pmra"], t["pmdec"], t["C"]), PDFs) # THIS TAKES VERY VERY LONG
g.add_column(t, RVs, "v_r")

t["pml"], t["pmb"] = zip(*[g.ICRS_to_galactic.transformProperMotions(p,th,mp,mth) for p,th,mp,mth in zip(t["ra_rad"], t["dec_rad"], t["pmra"], t["pmdec"])])
x,y,z,U,V,W = g.tophase(t["l"], t["b"], t["parallax"], t["pml"], t["pmb"], t["r_v"])
x.name, y.name, z.name, U.name, V.name, W.name = "x", "y", "z", "U", "V", "W"
t.add_columns((x,y,z,U,V,W))
print "Calculated individual stars' velocities"

#g.io.write_table_without_arrays(t, saveto="tgas/2MASS/TGAS_2MASS_afterXD.dat", exclude_cols=["solution_id", "random_index", "astrometric_n_obs_al", "astrometric_n_obs_ac", "astrometric_n_good_obs_al", "astrometric_n_good_obs_ac", "astrometric_n_bad_obs_al", "astrometric_n_bad_obs_ac", "astrometric_delta_q", "astrometric_excess_noise", "astrometric_excess_noise_sig", "astrometric_primary_flag", "astrometric_relegation_factor", "astrometric_weight_al", "astrometric_weight_ac", "astrometric_priors_used", "matched_observations", "duplicated_source", "scan_direction_strength_k1", "scan_direction_strength_k2", "scan_direction_strength_k3", "scan_direction_strength_k4", "scan_direction_mean_k1", "scan_direction_mean_k2", "scan_direction_mean_k3", "scan_direction_mean_k4", "phot_g_n_obs", "phot_g_mean_flux", "phot_g_mean_flux_error"])
g.io.write_table_with_separate_arrays(t, saveto_folder="tgas/2MASS/results/", exclude_cols=["solution_id", "random_index", "astrometric_n_obs_al", "astrometric_n_obs_ac", "astrometric_n_good_obs_al", "astrometric_n_good_obs_ac", "astrometric_n_bad_obs_al", "astrometric_n_bad_obs_ac", "astrometric_delta_q", "astrometric_excess_noise", "astrometric_excess_noise_sig", "astrometric_primary_flag", "astrometric_relegation_factor", "astrometric_weight_al", "astrometric_weight_ac", "astrometric_priors_used", "matched_observations", "duplicated_source", "scan_direction_strength_k1", "scan_direction_strength_k2", "scan_direction_strength_k3", "scan_direction_strength_k4", "scan_direction_mean_k1", "scan_direction_mean_k2", "scan_direction_mean_k3", "scan_direction_mean_k4", "phot_g_n_obs", "phot_g_mean_flux", "phot_g_mean_flux_error"])
print "Written file"

g.gplot.density(t["U"], t["V"], title = "TGAS-2MASS cross section $U$ vs. $V$", xlabel="$U$ (km/s)", ylabel="$V$ (km/s)", saveto="tgas/2MASS/TGAS_2MASS_XD_U_V.png", r=((-130,120),(-120,60)), bins = 200)
g.gplot.density(t["U"], t["W"], title = "TGAS-2MASS cross section $U$ vs. $W$", xlabel="$U$ (km/s)", ylabel="$W$ (km/s)", saveto="tgas/2MASS/TGAS_2MASS_XD_U_W.png", r=((-130,120),(-70,70)), bins = 200)
g.gplot.density(t["V"], t["W"], title = "TGAS-2MASS cross section $V$ vs. $W$", xlabel="$V$ (km/s)", ylabel="$W$ (km/s)", saveto="tgas/2MASS/TGAS_2MASS_XD_V_W.png", r=((-130,120),(-70,70)), bins = 200)

g.gplot.density(t["x"], t["y"], title = "TGAS-2MASS cross section $x$ vs. $y$", xlabel="$x$ (pc)", ylabel="$y$ (pc)", saveto="tgas/2MASS/TGAS_2MASS_x_y.png", bins = 200)
g.gplot.density(t["x"], t["z"], title = "TGAS-2MASS cross section $x$ vs. $z$", xlabel="$x$ (pc)", ylabel="$z$ (pc)", saveto="tgas/2MASS/TGAS_2MASS_x_z.png", bins = 200)
g.gplot.density(t["y"], t["z"], title = "TGAS-2MASS cross section $y$ vs. $z$", xlabel="$y$ (pc)", ylabel="$z$ (pc)", saveto="tgas/2MASS/TGAS_2MASS_y_z.png", bins = 200)
print "Plotted individual stars' velocities and positions"

B = (1000,0.02); range_binover=(-0.5, 3.7)
(meansU, sigmasU), binsU, highestU, sigmaBV = g.mapping.map_over_bins((g.np.mean, g.mean_sigma), t["g_min_ks"], B, -t["U"], range_binover=range_binover)
(meansV, sigmasV), binsV, highestV, sigmaBV = g.mapping.map_over_bins((g.np.mean, g.mean_sigma), t["g_min_ks"], B, -t["V"], range_binover=range_binover)
(meansW, sigmasW), binsW, highestW, sigmaBV = g.mapping.map_over_bins((g.np.mean, g.mean_sigma), t["g_min_ks"], B, -t["W"], range_binover=range_binover)
(S2, sigmaS2), binsS2, highestS2, sigmaBV = g.mapping.map_over_bins((g.S2, g.sigma_S2), t["g_min_ks"], B, t["ra_rad", "dec_rad", "U", "V", "W"], range_binover=range_binover)
bins = [(binsV, binsS2), binsU, binsW]
means = [(meansV, g.np.sqrt(S2)), meansU, meansW]
highests = [(highestV, highestS2), highestU, highestW]
sigmas = [(sigmasV, None), sigmasU, sigmasW]
kwargsU = {"ylabel": "$U$ [km/s]", "ylim": (-5,5)}
kwargsV = {"ylabel": "$V, S$ [km/s]", "ylim": (-5,48)}
kwargsW = {"ylabel": "$W$ [km/s]", "ylim": (-5,5), "xlabel": "$G - K$ [mag]", "xlim": range_binover}
kwargss = [kwargsV, kwargsU, kwargsW]
g.gplot.binned_1d_multi(bins, means, highests, sigmas_list = sigmas, kwargs_list = kwargss, sharex=True, gridspec_kw={"height_ratios": (2,1,1)}, saveto="tgas/2MASS/tgas_2mass_UVW_colour.png")

Umean = -t["U"].mean() ; Umeanerr = g.mean_sigma(t["U"])
Wmean = -t["W"].mean() ; Wmeanerr = g.mean_sigma(t["W"])
Vcoe, Vcov = g.np.polyfit(S2, meansV, 1, cov=True)
Vat0 = Vcoe[1] ; Vat0err = g.np.sqrt(Vcov[1,1])
print "Solar velocity:"
print "U: {0:.2f} +- {1:.2f}".format(Umean, Umeanerr)
print "V: {0:.2f} +- {1:.2f}".format(Vat0, Vat0err)
print "W: {0:.2f} +- {1:.2f}".format(Wmean, Wmeanerr)

g.gplot.S2UVW(S2, [meansU,meansV,meansW], [sigmasU,sigmasV,sigmasW], lines=[(0, Umean), Vcoe, (0, Wmean)], saveto="tgas/2MASS/S2UVW.png")

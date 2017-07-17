"""
Olivier Burggraaff
"""

import gaia_fc as g
print "Now loading Hipparcos table"
try:
    assert t
except:
    t = g.hip.read_multiple_data_files_as_tables("hip/hip_bovy2.vot")
    g.hip.add_C(t)
    print "Table loaded"

f = "hip/results/"
amps_xd = g.np.load(f+"amplitudes.npy")
means_xd = g.np.load(f+"means.npy")
covs_xd = g.np.load(f+"covariances.npy")

PDFs = map(g.pdf.multivariate, means_xd, covs_xd, amps_xd)

f, axs = g.gplot.plt.subplots(1,3,sharey=True, figsize=(20,7), tight_layout = True)
for H, vr, evr, ax in zip((50568, 10810, 87552), (-5.9,2.9,-19.7), (0.3,0.2,0.2), axs):
    print H,
    try:
        row = t[t["hip"] == H][0]
        print ""
    except IndexError:
        print "-----"
    RV = g.radial_velocity_distribution(PDFs, row["ra_rad"], row["de_rad"], row["plx"], row["pm_ra"], row["pm_de"], row["C"], x = g.np.arange(-150, 150, 0.1), v_r_err = evr)
    best_c = RV[0][RV[1].argmax()] ; H_c = g.pdf.entropy(RV[0], RV[1])
    best_m = RV[0][RV[2].argmax()] ; H_m = g.pdf.entropy(RV[0], RV[2])
    cdf = g.np.cumsum(RV[1:], axis=1) / 10.
    lower_c = RV[0][g.np.where(cdf[0] >= 0.025)[0][0]]
    lower_m = RV[0][g.np.where(cdf[1] >= 0.025)[0][0]]
    upper_c = RV[0][g.np.where(cdf[0] < 0.975)[0][-1]]
    upper_m = RV[0][g.np.where(cdf[1] < 0.975)[0][-1]]
    g.gplot.plotvrax(ax, RV[0], RV[1:], name = "HIP {0}".format(H), HRV = vr, eHRV = evr, lower_c = lower_c, upper_c = upper_c, H_c = H_c, H_m = H_m, xlabel="$v_r$ (km s$^{-1}$)", xlim=(-80,80))

axs[0].set_ylabel("$p(v_r)$")
for ax in axs:
    ax.set_xticks((-50, 0, 50))
f.savefig("hip_vr.png")

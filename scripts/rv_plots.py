from matplotlib import rcParams
import gaia_fc as g
rcParams.update({"font.size": "12"})
np = g.np
g.np.random.seed(1020202)
summary = g.io.read("tgas/summary_rave.txt", format="fixed_width")
summary.sort("ID")

indices = g.np.random.randint(0, len(summary), size = 6)
rows = summary[indices]

x = np.load("tgas/vr/xaxis.npy")

pm = 0
f, axs = g.gplot.plt.subplots(2, 3, sharey=True, sharex=True, figsize=(10,7), tight_layout = True)
for row, ax in zip(rows, axs.ravel()):
    print row["ID"],
    RV = np.load("tgas/vr/{0}.npy".format(row["ID"]))
    pm = max(pm, RV.max())
    g.gplot.plotvrax(ax, x, RV, HRV = row["HRV"], eHRV = row["eHRV"], lower_c = row["lower_c"], upper_c = row["upper_c"], H_c = row["H_c"], H_m = row["H_m"], xlim=(-80,80))
    
for ax in axs[1]:
    ax.set_xlabel("$v_r$ (km s$^{-1}$)")
for ax in axs[:,0]:
    ax.set_ylabel("$p(v_r)$")
    ax.set_ylim(0, 1.05*pm)
axs[0,0].set_yticks(ax.get_yticks()[::2])
axs[0,0].set_xticks([-50,0,50])
f.savefig("gaia_vr_example.png")
g.gplot.plt.close()

diff = g.np.abs(summary["HRV"] - summary["best_c"])
sigmas = (summary["upper_c"] - summary["lower_c"])/4.
ds = diff/sigmas

f, axs = g.gplot.plt.subplots(1,2, figsize=(10,5), tight_layout=True)
axs[0].hist(diff, bins = g.np.arange(0, 200, 5), normed=True)
axs[0].set_xlim(0, 75)
axs[0].set_xlabel("$|v_{r,\,RAVE} - v_{r,\,Gaia}|$ (km s$^{-1}$)")
axs[0].set_ylabel("Frequency")
axs[0].set_title("Mean: {mean:.0f} ; Median: {median:.0f} km".format(mean = g.np.mean(diff), median = g.np.median(diff))+" s$^{-1}$")
print "Fraction of stars not in absolute plot: {0:.1f}".format(100. * g.np.where(diff > 75)[0].shape[0] / float(len(diff)))
axs[1].hist(ds, bins = g.np.arange(0, 10, 0.25), normed=True)
axs[1].set_xlim(0, 5)
axs[1].set_xlabel("$|v_{r,\,RAVE} - v_{r,\,Gaia}| / \sigma_{Gaia}$")
axs[1].yaxis.tick_right()
axs[1].set_title("Mean: {mean:.1f} ; Median: {median:.1f} $\sigma$".format(mean = g.np.mean(ds), median = g.np.median(ds)))
print "Fraction of stars not in sigma plot: {0:.1f}".format(100. * g.np.where(ds > 5)[0].shape[0] / float(len(ds)))
f.savefig("rave_vs_2mass_2.png")

"""
g.gplot.plt.figure(figsize=(5,5))
g.gplot.plt.hist(diff, bins = g.np.arange(0, 200, 5), normed=True)
g.gplot.plt.xlim(0, 75)
g.gplot.plt.xlabel("$|v_{r,\,RAVE} - v_{r,\,Gaia}|$ (km s$^{-1}$)")
g.gplot.plt.ylabel("Frequency")
g.gplot.plt.title("Difference between observed and predicted $v_r$\nMean: {mean:.0f} ; Median: {median:.0f} km".format(mean = g.np.mean(diff), median = g.np.median(diff))+"s$^{-1}$")
print "Fraction of stars not in this plot: {0:.1f}".format(100. * g.np.where(diff > 75)[0].shape[0] / float(len(diff)))
g.gplot.plt.savefig("rave_vs_2mass.png")
g.gplot.plt.close()


g.gplot.plt.figure(figsize=(5,5))
g.gplot.plt.hist(ds, bins = g.np.arange(0, 10, 0.25), normed=True)
g.gplot.plt.xlim(0, 5)
g.gplot.plt.xlabel("$|v_{r,\,RAVE} - v_{r,\,Gaia}| / \sigma_{Gaia}$")
g.gplot.plt.ylabel("Frequency")
g.gplot.plt.title("Difference between observed and predicted $v_r$\nMean: {mean:.1f} ; Median: {median:.1f} $\sigma$".format(mean = g.np.mean(ds), median = g.np.median(ds)))
print "Fraction of stars not in this plot: {0:.1f}".format(100. * g.np.where(ds > 5)[0].shape[0] / float(len(ds)))
g.gplot.plt.savefig("rave_vs_2mass_R.png")
g.gplot.plt.close()
"""

print "Median 95% width:", g.np.median(sigmas*4.)

print "Within 95%:", 100.* g.np.where((summary["HRV"] > summary["lower_c"]) & (summary["HRV"] < summary["upper_c"]))[0].shape[0] / float(len(summary))

print "Mean error on HRV in RAVE:", g.np.mean(summary["eHRV"])

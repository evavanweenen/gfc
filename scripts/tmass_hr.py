import gaia_fc as g
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams.update({"font.size": 25, "figure.figsize": (10, 10)})

t = g.io.read_csv("tgas/2MASS/plx_small_noMS.csv")
print " --- read withoutMS" 
g.gplot.density(t["g_min_ks"], t["g_mag_abs"], saveto="tgasHRwithoutMS.png", xlabel="$G - K$", ylabel="$G$", flip="y", xlim=(-0.5, 4.0), ylim=(11, -1.8), cb=True, bins = 175, xticks=range(0, 5), title = "Final TGAS sample", histkwargs = {"vmin":1, "vmax":1200})
print "plotted"

t = g.io.read_votable("tgas/2MASS/plx_small.vot")
print " --- read withMS"
g.gplot.density(t["g_min_ks"], t["g_mag_abs"], saveto="tgasHRwithMS.png", xlabel="$G - K$", ylabel="$G$", flip="y", xlim=(-0.5, 4.0), ylim=(11, -1.8), cb=True, bins = 175, xticks=range(0, 5), title = "TGAS sample before main sequence cut", histkwargs = {"vmin":1, "vmax":1200})
print "plotted"

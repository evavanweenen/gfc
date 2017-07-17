"""
Find the optimal values of K, w for several XD models using two tests:
- XD internal likelihood
- predicted/observed radial velocities
"""

import gaia_fc as g
import matplotlib
matplotlib.rcParams.update({'font.size': 26, "figure.figsize": (10, 5)})
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator

sqwrange = g.np.arange(0.5, 3.75, 0.25)
wrange = sqwrange**2.
Krange = range(4, 21)

Ltable = g.np.tile(g.np.nan, (len(Krange), len(wrange)))
LtableR = g.np.tile(g.np.nan, (len(Krange), len(wrange)))

for i,w in enumerate(wrange):
    for j,K in enumerate(Krange):
        f = "tgas/2MASS/XD_K_w/w_{w}/K_{K}/L".format(w = w, K = K)
        try:
            L = g.np.loadtxt(f)
            Ltable[j,i] = L
        except:
            pass
        try:
            L_rave = g.np.loadtxt(f+"_rave")
            if L_rave < -1000:
                LtableR[j,i] = L_rave
        except:
            pass

Ltable += 27
LtableR += 247000

bestL = g.np.unravel_index(g.np.nanargmax(Ltable), Ltable.shape)
bestLR = g.np.unravel_index(g.np.nanargmax(LtableR), LtableR.shape)
bestKL = Krange[bestL[0]] ; bestwL = wrange[bestL[1]]
print "L: best K = {K} ; best w = {w}".format(K = bestKL, w = bestwL)
bestKLR = Krange[bestLR[0]] ; bestwLR = wrange[bestLR[1]]
print "LR: best K = {K} ; best w = {w}".format(K = bestKLR, w = bestwLR)

f, axs = g.gplot.plt.subplots(1,2,sharey=True,tight_layout=True, figsize=(15,7))
im0 = axs[0].imshow(Ltable.T, cmap=g.gplot.plt.cm.viridis, extent=(Krange[0]-0.5, Krange[-1]+0.5, sqwrange[0]-0.125, sqwrange[-1]+0.125), aspect='auto', interpolation='none', origin="lower", vmin=-0.273)
axs[0].scatter(bestKL, g.np.sqrt(bestwL), color="black", marker="o", s=250)
im1 = axs[1].imshow(LtableR.T, cmap=g.gplot.plt.cm.viridis, extent=(Krange[0]-0.5, Krange[-1]+0.5, sqwrange[0]-0.125, sqwrange[-1]+0.125), aspect='auto', interpolation='none', origin="lower", vmin=125)
axs[1].scatter(bestKLR, g.np.sqrt(bestwLR), color="black", marker="o", s=250)
axs[0].set_xlabel("$K$") ; axs[1].set_xlabel("$K$")
axs[0].set_ylabel("$\sqrt{w}$ (km s$^{-1}$)")
axs[0].set_xlim(Krange[0]-0.5, Krange[-1]+0.5)
axs[1].set_xlim(Krange[0]-0.5, Krange[-1]+0.5)
axs[0].set_ylim(sqwrange[0]-0.125, sqwrange[-1]+0.125)

div0 = make_axes_locatable(axs[0])
div1 = make_axes_locatable(axs[1])
cax0 = div0.append_axes("top", size="8%", pad=0.05)
cbar0 = f.colorbar(im0, cax=cax0, format="%.3f", orientation="horizontal")
cax0.xaxis.tick_top()
cbar0.set_ticks(g.np.arange(-0.273, -0.265, 0.002))
cax1 = div1.append_axes("top", size="8%", pad=0.05)
cbar1 = f.colorbar(im1, cax=cax1, format="%.0f", orientation="horizontal")
cax1.xaxis.tick_top()
cbar1.set_ticks(g.np.arange(150, 650, 100))

f.savefig("L.png")


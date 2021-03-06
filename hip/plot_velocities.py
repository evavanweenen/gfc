"""
Plot individual components and total velocity distribution of XD model
"""

import gaia_fc as g
from sys import argv
from matplotlib import pyplot as plt
from matplotlib import gridspec as gs
from matplotlib.patches import Ellipse
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

from scipy.stats import chi2

import numpy as np

#K = int(argv[1]); w = float(argv[2])
K = 19 ; w = 9.0

volume = 0.6827

def text(ax1):
    ax1.text(60, -90, "Arcturus")
    ax1.plot([8, 57.5], [-105, -90], c='k', lw=2)
    ax1.text(-80, -110, "Halo")
    ax1.text(50, 40, "Sirius/UMa")
    ax1.plot([49, 5], [40, 5], c='k', lw=2)
    ax1.text(-100, 45, "Coma Berenices")
    ax1.plot([-70, -10], [42, -5], c='k', lw=2)
    ax1.text(-120, 34, "NGC 1901")
    ax1.plot([-100, -25], [31, -12], c='k', lw=2)
    ax1.text(-120, 0, "Hyades")
    ax1.plot([-110, -45], [-3, -17], c='k', lw=2)
    ax1.text(90, -50, "Pleiades")
    ax1.plot([87, -15], [-45, -20], c='k', lw=2)
    ax1.text(-125, -42, "Hercules")
    ax1.plot([-93.5, -28], [-40, -42], c='k', lw=2)
    
def eigsorted(cov):
    """
    http://www.nhsilbert.net/source/2014/06/bivariate-normal-ellipse-plotting-in-python/
    """
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]
def draw_PDF_ellipse(ax, amp, mean, cov, xyz, volume = 0.6827, **kwargs):
    assert xyz in ("xy", "xz", "yz"), "gaia_fc.gplot.draw_PDF_ellipse: parameter `xyz` must be one of ('xy', 'xz', 'yz'); received value {0}".format(xyz)
    choose_from = {"xy": [0, 1], "xz": [0, 2], "yz": [1, 2]}
    pair_int = choose_from[xyz]
    mean_ = mean[pair_int]
    cov_ = cov[pair_int][:, pair_int]

    vals, vecs = eigsorted(cov_)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * np.sqrt(chi2.ppf(volume,2)) * np.sqrt(vals)
    ell = Ellipse(xy = mean_, width = width, height = height, angle = theta, facecolor = "None", linewidth = 1 + 7.5*amp, **kwargs)
    ax.add_artist(ell)

def set_axes(ax1, ax2, ax3):
    ax1.set_xlim(-130, 120) ; ax1.set_ylim(-120, 60)
    ax2.set_xlim(-130, 120) ; ax2.set_ylim(-70, 70)
    ax3.set_xlim(-130, 120) ; ax3.set_ylim(-70, 70)
    ax1.set_xlabel("$U$ (km s$^{-1}$)") ; ax1.xaxis.set_label_position("top")
    ax1.set_ylabel("$V$ (km s$^{-1}$)") ; ax1.yaxis.set_label_position("left")
    ax2.set_xlabel("$U$ (km s$^{-1}$)") ; ax2.xaxis.set_label_position("bottom")
    ax2.set_ylabel("$W$ (km s$^{-1}$)") ; ax2.yaxis.set_label_position("left")
    ax3.set_xlabel("$V$ (km s$^{-1}$)") ; ax3.xaxis.set_label_position("bottom")
    ax3.set_ylabel("$W$ (km s$^{-1}$)") ; ax3.yaxis.set_label_position("right")
    ax1.xaxis.tick_top()
    ax3.yaxis.tick_right()
    ax1.set_xticks((-100, -50, 0, 50, 100))
    ax1.set_yticks((-100, -50, 0, 50))

print "Doing Hipparcos"
amps =  np.load(r"hip/results/amplitudes.npy")
means = np.load(r"hip/results/means.npy")
covs = np.load(r"hip/results/covariances.npy")
evalxyz = np.load("hip/results/eval.npy")

f = plt.figure(figsize=(8,8), tight_layout = True)

plt.suptitle("Gaussian components of velocity distribution fitted to Hipparcos data")
gs1 = gs.GridSpec(3, 4)
gs1.update(left=0.1, right=0.9, wspace=0.05)
ax1 = plt.subplot(gs1[:-1, :])
ax2 = plt.subplot(gs1[-1, :2])
ax3 = plt.subplot(gs1[-1, 2:])
set_axes(ax1, ax2, ax3)

for a, m, c in zip(amps, means, covs):
    draw_PDF_ellipse(ax1, a, m, c, "xy", edgecolor="0.5")
    draw_PDF_ellipse(ax2, a, m, c, "xz", edgecolor="0.5")
    draw_PDF_ellipse(ax3, a, m, c, "yz", edgecolor="0.5")

text(ax1)

f.savefig("PDFs_hip.png")
plt.close(f)

print " - Components plotted"

levels = np.logspace(-6, -2.7, 12)

evalxy = evalxyz.sum(2)
evalxz = evalxyz.sum(1)
evalyz = evalxyz.sum(0)
f = plt.figure(figsize=(8,8), tight_layout = True)

plt.suptitle("Velocity distribution fitted to Hipparcos data")
gs1 = gs.GridSpec(3, 4)
gs1.update(left=0.1, right=0.9, wspace=0.05)
ax1 = plt.subplot(gs1[:-1, :])
ax2 = plt.subplot(gs1[-1, :2])
ax3 = plt.subplot(gs1[-1, 2:])
ax1.contour(evalxy.T, extent = [-140, 140, -130, 130], levels = levels, colors = '0.5')
ax2.contour(evalxz.T, extent = [-140, 140, -72, 72], levels = levels, colors = '0.5')
ax3.contour(evalyz.T, extent = [-130, 130, -72, 72], levels = levels, colors = '0.5')
set_axes(ax1, ax2, ax3)
text(ax1)
f.savefig("Eval_hip.png")
plt.close(f)

print " - Total plotted"

s = []
s.append(r"\begin{table}[h]")
s.append(r"\centering")
s.append(r"\caption{Parameters of the components of the model with $K = 10$ and $w = 4$ km$^2$ s$^{-2}$}")
s.append(r"\label{t:hc}")
s.append(r"\begin{adjustbox}{center}")
s.append(r"\begin{tabular}{rlrrrrrrrrr}")
s.append(r"\hline\hline")
s.append(r"$j$ & $\alpha$ & $\unit x \tran \vec m$ & $\unit y \tran \vec m$ & $\unit z \tran \vec m$ & $\unit x \tran \mat V \unit x$ & $\unit y \tran \mat V \unit y$ & $\unit z \tran \mat V \unit z$ & $\unit x \tran \mat V \unit y$  & $\unit x \tran \mat V \unit z$  & $\unit y \tran \mat V \unit z$ \\")
s.append(r" & & (km s$^{-1}$)& (km s$^{-1}$)& (km s$^{-1}$)& (km$^2$ s$^{-2}$)& (km$^2$ s$^{-2}$)& (km$^2$ s$^{-2}$)& (km$^2$ s$^{-2}$)& (km$^2$ s$^{-2}$)& (km$^2$ s$^{-2}$) \\")
s.append(r"\hline")
for j, a, m, c in zip(range(len(amps)), amps, means, covs):
    l = "{j} & {a:.4f} & {xv:.2f} & {yv:.2f} & {zv:.2f} & {xvx:.2f} & {yvy:.2f} & {zvz:.2f} & {xvy:.2f} & {xvz:.2f} & {yvz:.2f}"
    l = l.format(j = j+1, a = a, xv = m[0], yv = m[1], zv = m[2], xvx = c[0,0], yvy = c[1,1], zvz = c[2,2], xvy = c[0,1], xvz = c[0,2], yvz = c[1,2])
    s.append(l + r"\\")

s.append(r"\end{tabular}")
s.append(r"\end{adjustbox}")
s.append(r"\end{table}")

s2 = reduce(lambda x, y: x + "\n" + y, s)
print >> open("hip_comp.tex", 'w'), s2

print " - Table written"

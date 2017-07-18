"""
Plot individual components and total velocity distribution of XD model

TODO: Move to gfc.gplot
"""

import gfc
from matplotlib import pyplot as plt
from matplotlib import gridspec as gs
from matplotlib.patches import Ellipse
import matplotlib as mpl

from scipy.stats import chi2

import numpy as np

from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("xd_results_folder", help = "Folder that contains the XD results")
parser.add_argument("saveto_folder", help = "Folder in which plot will be saved")
parser.add_argument("-n", "--nrlevels", help = "Number of levels in total plot", type = int, default = 12)
parser.add_argument("--minlevel", help = "(Log) minimum level in total plot", type = float, default = -6.0)
parser.add_argument("--maxlevel", help = "(Log) maximum level in total plot", type = float, default = -2.7)
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

volume = 0.6827
levels = np.logspace(args.minlevel, args.maxlevel, args.nrlevels)

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

def eigsorted(cov):
    """
    http://www.nhsilbert.net/source/2014/06/bivariate-normal-ellipse-plotting-in-python/
    """
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]
def draw_PDF_ellipse(ax, amp, mean, cov, xyz, volume = 0.6827, **kwargs):
    choose_from = {"xy": [0, 1], "xz": [0, 2], "yz": [1, 2]}
    pair_int = choose_from[xyz]
    mean_ = mean[pair_int]
    cov_ = cov[pair_int][:, pair_int]

    vals, vecs = eigsorted(cov_)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * np.sqrt(chi2.ppf(volume,2)) * np.sqrt(vals)
    ell = Ellipse(xy = mean_, width = width, height = height, angle = theta, facecolor = "None", linewidth = 1 + 7.5*amp, **kwargs)
    ax.add_artist(ell)

amps_xd = np.load(args.xd_results_folder+"amplitudes.npy")
means_xd = np.load(args.xd_results_folder+"means.npy")
covs_xd = np.load(args.xd_results_folder+"covariances.npy")

f = plt.figure(figsize=(8,8), tight_layout = True)

plt.suptitle("Gaussian components of velocity distribution fitted to Gaia data")
gs1 = gs.GridSpec(3, 4)
gs1.update(left=0.1, right=0.9, wspace=0.05)
ax1 = plt.subplot(gs1[:-1, :])
ax2 = plt.subplot(gs1[-1, :2])
ax3 = plt.subplot(gs1[-1, 2:])
set_axes(ax1, ax2, ax3)

for a, m, c in zip(amps_xd, means_xd, covs_xd):
    draw_PDF_ellipse(ax1, a, m, c, "xy", edgecolor="0.4")
    draw_PDF_ellipse(ax2, a, m, c, "xz", edgecolor="0.4")
    draw_PDF_ellipse(ax3, a, m, c, "yz", edgecolor="0.4")

text(ax1)

f.savefig(args.saveto_folder+"/PDFs.png")
plt.close(f)

PDFs = map(gfc.pdf.multivariate, means_xd, covs_xd, amps_xd)
evalxyz = gfc.pdf.eval_total_PDF(PDFs, [(-140,140), (-130,130), (-72,72)])
evalxy = evalxyz.sum(2)
evalxz = evalxyz.sum(1)
evalyz = evalxyz.sum(0)
f = plt.figure(figsize=(8,8), tight_layout = True)

plt.suptitle("Velocity distribution fitted to Gaia data")
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
f.savefig(args.saveto_folder+"/total.png")
plt.close(f)

print " - Total plotted"

s = []
s.append(r"\begin{table}[h]")
s.append(r"\centering")
s.append(r"\caption{Parameters of the components of the model")# with $K = "+str(K)+r"$ and $w = "+str(w)+"$ km$^2$ s$^{-2}$}")
s.append(r"\label{t:gc}")
s.append(r"\begin{adjustbox}{center}")
s.append(r"\begin{tabular}{rlrrrrrrrrr}")
s.append(r"\hline\hline")
s.append(r"$j$ & $\alpha$ & $\unit x \tran \vec m$ & $\unit y \tran \vec m$ & $\unit z \tran \vec m$ & $\unit x \tran \mat V \unit x$ & $\unit y \tran \mat V \unit y$ & $\unit z \tran \mat V \unit z$ & $\unit x \tran \mat V \unit y$  & $\unit x \tran \mat V \unit z$  & $\unit y \tran \mat V \unit z$ \\")
s.append(r" & & (km s$^{-1}$)& (km s$^{-1}$)& (km s$^{-1}$)& (km$^2$ s$^{-2}$)& (km$^2$ s$^{-2}$)& (km$^2$ s$^{-2}$)& (km$^2$ s$^{-2}$)& (km$^2$ s$^{-2}$)& (km$^2$ s$^{-2}$) \\")
s.append(r"\hline")
for j, a, m, c in zip(range(len(amps_xd)), amps_xd, means_xd, covs_xd):
    l = "{j} & {a:.5f} & {xv:.2f} & {yv:.2f} & {zv:.2f} & {xvx:.2f} & {yvy:.2f} & {zvz:.2f} & {xvy:.2f} & {xvz:.2f} & {yvz:.2f}"
    l = l.format(j = j+1, a = a, xv = m[0], yv = m[1], zv = m[2], xvx = c[0,0], yvy = c[1,1], zvz = c[2,2], xvy = c[0,1], xvz = c[0,2], yvz = c[1,2])
    s.append(l + r"\\")

s.append(r"\end{tabular}")
s.append(r"\end{adjustbox}")
s.append(r"\end{table}")

s2 = reduce(lambda x, y: x + "\n" + y, s)
print >> open(args.saveto_folder+"gaia_comp.tex", 'w'), s2

print " - Table written"

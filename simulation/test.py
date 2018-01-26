import gfc
from gfc import *

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LogNorm
from matplotlib import rc
rc('font',**{'family':'serif', 'serif':['Times']})
rc('text', usetex=True)
import copy as cp
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import Angle

import pygaia.astrometry.vectorastrometry as pg

import scipy as sp

parser = ArgumentParser()
parser.add_argument("--i", help = "Number of components", default = 1, type= int)
parser.add_argument("--L", help = "Length of data", default=1000, type = int)
args = parser.parse_args()

wmin = 0.5
wmax=5
wstep=.25
Kmin=1
Kmax=10
Kstep=1

args.L = 1000

init_amps = np.empty([10,10])
for i in xrange(10):
    amp = np.tile([1./(i+1)], i+1)
    rest = np.tile(0, 9-i)
    init_amps[i] = np.append(amp, rest)
    print "init_amps", init_amps

init_means = np.array([[10.], [20.], [30.], [40.], [50.], [60.], [70.], [80.], [90.], [100.]])
init_covs = np.array([[[0.1]], [[0.1]], [[0.1]], [[0.1]], [[0.1]], [[0.1]], [[0.1]], [[0.1]], [[0.1]], [[0.1]]])
print "init_means", init_means
print "init_covs", init_covs

comp = np.random.choice(args.i, args.L, p=init_amps[args.i-1][:args.i])
print "comp", comp

d = np.empty([args.L,1])
dcov = np.empty([args.L,1])
for l in xrange(args.L):
    d[l,0] = np.random.multivariate_normal(init_means[comp[l]], init_covs[comp[l]])
    dcov[l,0] = init_covs[comp[l]]

print "d", d
print "dcov", dcov

proj = np.tile(np.identity(1), (args.L,1))

plt.title("Histogram of K = {0}".format(args.i))
plt.hist(d, facecolor='green', normed=True, histtype='stepfilled', alpha=.3)
plt.savefig("../results/test/test/Histogram_K{K}_L{L}.png".format(K = args.i, L = args.L))
plt.show()

wrange = gfc.np.arange(wmin, wmax + wstep, wstep)**2.
Krange = range(Kmin, Kmax + Kstep, Kstep)

LLH_table = np.tile(np.nan, (len(Krange), len(wrange)))

for i,k in enumerate(Krange):
    for j,w in enumerate(wrange):
        input_amps = cp.copy(init_amps)
        input_means = cp.copy(init_means)
        input_covs = cp.copy(init_covs)
        print "init_amps = {0} ; init_means = {1} ; init_covs = {2}".format(input_amps[k-1][:k], input_means[:k], input_covs[:k])
        amps_xd, means_xd, covs_xd, LLH = gfc.XD(d, dcov, input_amps[k-1][:k], input_means[:k], input_covs[:k], projection=proj, w=w) #extreme deconvolution
        LLH_table[i,j] = LLH
        print "K = {0} ; w = {1} ; LLH = {2}".format(k, w, LLH)
        print "amps_xd = {0} ; means_xd = {1} ; covs_xd = {2}".format(amps_xd, means_xd, covs_xd)


print "LLH_table", LLH_table

best_LLH = np.unravel_index(np.nanargmax(LLH_table), LLH_table.shape)
print best_LLH
bestK = Krange[best_LLH[0]] ; bestw = wrange[best_LLH[1]]
print "L: best K = {K} ; best w = {w}".format(K = bestK, w = bestw)


xmin = Krange[0]-0.5
xmax = Krange[-1]+0.5
ymin = np.sqrt(wrange[0])-0.125
ymax = np.sqrt(wrange[-1])+0.125
xlabel = "$K$"
ylabel = "$\sqrt{w}$ (km s$^{-1}$)"
title = "Likelihood K={K}, L={L}, bestw={bestw}, bestK={bestK}".format(K = args.i, L = args.L, bestw = bestw, bestK = bestK)

fig = plt.figure(figsize=(5,4))
im = plt.imshow(LLH_table.T, cmap=plt.cm.viridis, extent=(xmin, xmax, ymin, ymax), aspect='auto', interpolation='none', origin="lower")
plt.scatter(bestK, np.sqrt(bestw), color="black", marker="o", s=250)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.title(title)
cbar = plt.colorbar(im, format="%.3f", orientation="vertical")

file_name = "../results/test/test/K{K}_L{L}_bestw{bestw}_bestK{bestK}.png".format(K = args.i, L = args.L, bestw = bestw, bestK=bestK)
plt.savefig(file_name)


"""
Fit a model to TGAS data using XD.
"""
import gfc
from gfc import *
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LogNorm
import pygaia.astrometry.vectorastrometry as pg

from matplotlib import rc
rc('font',**{'family':'serif', 'serif':['Times']})
rc('text', usetex=True)

from gfc import ArgumentParser
parser = ArgumentParser()

parser.add_argument("root_folder", help = "Folder in which results will be saved")
parser.add_argument("--K", help = "True value of K", type=int)
parser.add_argument("--L", help = "Number of data points", type=int)
parser.add_argument("--Kmin", help = "Lowest value of K to test", default = 1, type = int)
parser.add_argument("--Kmax", help = "Highest value of K to test", default = 10, type = int)
parser.add_argument("--Kstep", help = "Step between K values", default = 1, type = int)
parser.add_argument("--wmin", help = "Lowest value of SQRT(w) to test", default = 0.5, type = float)
parser.add_argument("--wmax", help = "Highest value of SQRT(w) to test", default = 5, type = float)
parser.add_argument("--wstep", help = "Step between SQRT(w) values", default = 0.25, type = float)
parser.add_argument("--wsigma2", help = "Sigma you want to use", default = 3., type = float)
parser.add_argument("--offset", help = "Offset of data", default = 27, type = int)
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

if args.verbose:
    print "Finished parsing arguments"

"""
Create w and K array
"""    
wrange = gfc.np.arange(args.wmin, args.wmax + args.wstep, args.wstep)**2.
Krange = range(args.Kmin, args.Kmax + args.Kstep, args.Kstep)

"""
Create separate folder for each different model (so for each w and K combination)   .
Olivier: Keep in mind that this is not the best solution for the w-K likelihood plot.
"""
save_folder = args.root_folder + "simulated_K{K}/L{L}".format(K = args.K, L = args.L)
print "save_folder:", save_folder

for w in wrange:
    f_w = "{f}/w_{w}".format(f = save_folder, w = w)
    if not gfc.isdir(f_w):
        gfc.mkdir(f_w)
    for K in Krange:
        f = "{f}/w_{w}/K_{K}".format(f = save_folder, w = w, K = K)
        if not gfc.isdir(f):
            gfc.mkdir(f)

if args.verbose:
    print "Finished creating folders"


"""
Read data
"""
data_file = args.root_folder + "simulated_K{K}/simulated_data_xv_K{K}_L{L}.npy".format(K = args.K, L = args.L)
print "data_file:", data_file
time_before_loading = gfc.time() #times how long programme is running
t = gfc.io.read_csv(data_file)
gfc.remove_unused_columns(t)
time_after_loading = gfc.time()

if args.verbose:
    print "Finished loading data in {0:.1f} seconds".format(time_after_loading - time_before_loading)


"""
Calculate vectors and covariance matrix
"""
"""
time_before_matrices = gfc.time()
gfc.add_rad(t)
gfc.matrix.add_w(t)
gfc.matrix.add_A(t)
gfc.matrix.add_R(t)
gfc.tgas.add_C(t)
gfc.matrix.add_Q(t)
gfc.matrix.add_S(t)
time_after_matrices = gfc.time()
if args.verbose:
    print "Added w, A, R, C, Q, S in {0:.1f} seconds".format(time_after_matrices - time_before_matrices)
"""

"""
XD cannot handle csv format, so convert to numpy tables
"""
#warr = gfc.XD_arr(t, "w1", "w2", "w3")
#wcov = gfc.XD_arr(t, "S") ; wcov[:, 0, 0] = 1e15 #hardcode: pretend first component (vrad) is 0, and make covariance very high so it is ignored
#proj = gfc.XD_arr(t, "R")

"""
err = np.empty([6])
err = pg.astrometryToPhaseSpace(0.,0.,.3,1.,1.,1e15)
wcov_one = np.array([[err[3],0,0],[0,err[4],0],[0,0,err[5]]])
"""
warr = gfc.XD_arr(t, "w1", "w2", "w3")
wcov_one = np.array([[args.wsigma2,0,0],[0,args.wsigma2,0],[0,0,args.wsigma2]])
wcov = np.tile(wcov_one, (args.L,1,1))
proj = np.tile(np.identity(3), (args.L,1,1))
print "warr", warr
print "wcov", wcov
print "proj", proj


"""
Perform XD and write to file
"""
LLH_table = np.tile(np.nan, (len(Krange), len(wrange)))
for i,k in enumerate(Krange):
    for j,w in enumerate(wrange):
        f = "{folder}/w_{w}/K_{k}".format(folder = save_folder, w = w, k = k)
        init_amps = gfc.io.load(args.root_folder + "/initial_values/initial_amps_K{k}.npy".format(k = k))
        init_means = gfc.io.load(args.root_folder + "/initial_values/initial_means_K{k}.npy".format(k = k))
        init_covs = gfc.io.load(args.root_folder + "/initial_values/initial_covs_K{k}.npy".format(k = k))
        print "Initial values before XD: initamps = {0} ; initmeans = {1} ; initcovs = {2}".format(init_amps, init_means, init_covs)
        amps_xd, means_xd, covs_xd, LLH = gfc.XD(warr, wcov, init_amps, init_means, init_covs, projection=proj, w=w) #extreme deconvolution
        LLH_table[i,j] = LLH
        #print "Initial values after XD: initamps = {0} ; initmeans = {1} ; initcovs = {2}".format(init_amps, init_means, init_covs)
        #print "K = {0} ; w = {1} ; LLH = {2}".format(k, w, LLH)
        print "LLH = {0} ; amps_xd = {1} ; means_xd = {2} ; covs_xd = {3}".format(LLH, amps_xd, means_xd, covs_xd)
        print >> open("{0}/L".format(f), 'w'), LLH
        gfc.io.save_PDFs(amps_xd, means_xd, covs_xd, f)

#LLH_table += args.offset

print LLH_table

"""
Find bestK and bestw
"""        
best_LLH = np.unravel_index(np.nanargmax(LLH_table), LLH_table.shape)
bestK = Krange[best_LLH[0]] ; bestw = wrange[best_LLH[1]]
print "L: best K = {K} ; best w = {w}".format(K = bestK, w = bestw)


"""
Plot bestK and bestw
"""
xmin = Krange[0]-0.5
xmax = Krange[-1]+0.5
ymin = np.sqrt(wrange[0])-0.125
ymax = np.sqrt(wrange[-1])+0.125
xlabel = "$K$"
ylabel = "$\sqrt{w}$ (km s$^{-1}$)"
title = "Likelihood K={K}, L={L}, bestw={bestw}, bestK={bestK}".format(K = args.K, L = args.L, bestw = bestw, bestK = bestK)

fig = plt.figure(figsize=(5,4))
im = plt.imshow(LLH_table.T, cmap=plt.cm.viridis, extent=(xmin, xmax, ymin, ymax), aspect='auto', interpolation='none', origin="lower")
plt.scatter(bestK, np.sqrt(bestw), color="black", marker="o", s=250)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.title(title)
cbar = plt.colorbar(im, format="%.3f", orientation="vertical")

file_name = "/L__xv_K{K}_L{L}_wmax{wmax}_bestw{bestw}_bestK{bestK}.png".format(K = args.K, L = args.L, wmax=args.wmax**2, bestw = bestw, bestK=bestK)
plt.savefig(save_folder + file_name)

"""
Fit a model to TGAS data using XD.
Note: several matrices must be pre-calculated, e.g. using the script `xd_calculate_matrices.py`
"""
import gfc

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("data_folder", help = "Folder that contains the data")
parser.add_argument("save_folder", help = "Folder in which results will be saved")
parser.add_argument("--init_amps", help = "File with initial estimates of amplitudes", default = "Bovy_parameters/Bovy_amps.npy")
parser.add_argument("--init_means", help = "File with initial estimates of means", default = "Bovy_parameters/Bovy_means.npy")
parser.add_argument("--init_covs", help = "File with initial estimates of covariances", default = "Bovy_parameters/Bovy_covs.npy")
parser.add_argument("--Kmin", help = "Lowest value of K to test", default = 4, type = int)
parser.add_argument("--Kmax", help = "Highest value of K to test", default = 20, type = int)
parser.add_argument("--Kstep", help = "Step between K values", default = 1, type = int)
parser.add_argument("--wmin", help = "Lowest value of SQRT(w) to test", default = 0.5, type = float)
parser.add_argument("--wmax", help = "Highest value of SQRT(w) to test", default = 3.5, type = float)
parser.add_argument("--wstep", help = "Step between SQRT(w) values", default = 0.25, type = float)
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

if args.verbose:
    print "Finished parsing arguments"

initial_amps = gfc.io.load(args.init_amps)
initial_means = gfc.io.load(args.init_means)
initial_covs = gfc.io.load(args.init_covs)

if args.verbose:
    print "Loaded initial estimates for Gaussian parameters"
    
wrange = gfc.np.arange(args.wmin, args.wmax + args.wstep, args.wstep)**2.
Krange = range(args.Kmin, args.Kmax + args.Kstep, args.Kstep)

for w in wrange:
    f_w = "{f}/w_{w}".format(f = args.save_folder, w = w)
    if not gfc.isdir(f_w):
        gfc.mkdir(f_w)
    for K in Krange:
        f = "{f}/w_{w}/K_{K}".format(f = args.save_folder, w = w, K = K)
        if not gfc.isdir(f):
            gfc.mkdir(f)

if args.verbose:
    print "Finished creating folders"

t = gfc.io.load_table_with_separate_arrays(saveto_folder = args.data_folder)
assert all(col in t.keys() for col in ["w1", "w2", "w3", "S", "R"])

if args.verbose:
    print "Finished loading data"

warr = gfc.XD_arr(t, "w1", "w2", "w3")
wcov = gfc.XD_arr(t, "S") ; wcov[:, 0, 0] = 1e15
proj = gfc.XD_arr(t, "R")

for K in Krange:
    for w in wrange:
        f = "{f}/w_{w}/K_{K}".format(f = args.save_folder, w = w, K = K)
        print "K = {0} ; w = {1}".format(K, w)
        amps_xd, means_xd, covs_xd, L = gfc.XD(warr, wcov, initial_amps[:K], initial_means[:K], initial_covs[:K], projection = proj, w = w)
        print >> open(f+"/L", 'w'), L
        gfc.io.save_PDFs(amps_xd, means_xd, covs_xd, f)

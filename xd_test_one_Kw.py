"""
Fit a model to TGAS data using XD.
Note: several matrices must be pre-calculated, e.g. using the script `xd_calculate_matrices.py`
"""
import gfc

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("data_folder", help = "Folder that contains the data")
parser.add_argument("save_folder", help = "Folder in which results will be saved")
parser.add_argument("-K", help = "K", default = 15, type = int)
parser.add_argument("-w", help = "w", default = 4.0, type = float)
parser.add_argument("--init_amps", help = "File with initial estimates of amplitudes", default = "Bovy_parameters/Bovy_amps.npy")
parser.add_argument("--init_means", help = "File with initial estimates of means", default = "Bovy_parameters/Bovy_means.npy")
parser.add_argument("--init_covs", help = "File with initial estimates of covariances", default = "Bovy_parameters/Bovy_covs.npy")
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

if args.verbose:
    print "Finished parsing arguments"

initial_amps = gfc.io.load(args.init_amps)
initial_means = gfc.io.load(args.init_means)
initial_covs = gfc.io.load(args.init_covs)
K = args.K
w = args.w

if args.verbose:
    print "Loaded initial estimates for Gaussian parameters"
    
t = gfc.io.load_table_with_separate_arrays(saveto_folder = args.data_folder)
assert all(col in t.keys() for col in ["w1", "w2", "w3", "S", "R"])

if args.verbose:
    print "Finished loading data"

warr = gfc.XD_arr(t, "w1", "w2", "w3")
wcov = gfc.XD_arr(t, "S") ; wcov[:, 0, 0] = 1e15
proj = gfc.XD_arr(t, "R")

print "K = {0} ; w = {1}".format(K, w)
amps_xd, means_xd, covs_xd, L = gfc.XD(warr, wcov, initial_amps[:K], initial_means[:K], initial_covs[:K], projection = proj, w = w)
print >> open("{0}/L".format(args.save_folder), 'w'), L
gfc.io.save_PDFs(amps_xd, means_xd, covs_xd, args.save_folder)

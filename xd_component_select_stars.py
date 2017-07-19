import gfc
import numpy as np

from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("data_file", help = "File that contains the stellar data")
parser.add_argument("xd_results_folder", help = "Folder that contains the XD results")
parser.add_argument("save_folder", help = "Folder in which to save results")
parser.add_argument("-v", "--verbose", action = "store_true")
parser.add_argument("-t", "--threshold", help = "-1 * Log-likelihood threshold for component membership", default = 9, type = int)
args = parser.parse_args()

t = gfc.io.read_csv(args.data_file)

if args.verbose:
    print "Finished loading data"

amps_xd = np.load(args.xd_results_folder+"amplitudes.npy")
means_xd = np.load(args.xd_results_folder+"means.npy")
covs_xd = np.load(args.xd_results_folder+"covariances.npy")

if args.verbose:
    print "Finished loading XD parameters"

gfc.tgas.add_rad(t)
gfc.tgas.add_w(t)
gfc.tgas.add_A(t)
gfc.tgas.add_R(t)
gfc.tgas.add_UVW(t)

if args.verbose:
    print "Calculated U, V, W"

Ls = gfc.pdf.loglikelihood_many_multiPDF(means_xd, covs_xd, t["UVW_vec"])

if args.verbose:
    print "Calculated likelihoods"

gfc.gplot.moving_group(t["ra"], t["dec"], t["l"], t["b"], t["U"], t["V"], t["W"], range(len(t)), [], bins = 100, saveto = "all_mg.png") 
for n, (a, m, c) in enumerate(zip(amps_xd, means_xd, covs_xd)):
    indices_in = np.where(Ls[:,n] > -args.threshold)[0]
    indices_ex = np.where(Ls[:,n] <= -args.threshold)[0]
    if args.verbose:
        print "{n:02d}: {l}".format(n=n, l=len(t_))
    gfc.gplot.moving_group(t_["ra"], t_["dec"], t_["l"], t_["b"], t_["U"], t_["V"], t_["W"], indices_in, indides_ex, \
        a = a, m = m, c = c, bins = 100, saveto = "{f}/{n:02d}_mg.png".format(f = args.save_folder, n = n+1)) 

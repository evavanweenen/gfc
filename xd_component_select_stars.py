import gfc
import numpy as np

from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("data_file", help = "File that contains the stellar data")
parser.add_argument("xd_results_folder", help = "Folder that contains the XD results")
parser.add_argument("save_folder", help = "Folder in which to save results")
parser.add_argument("-r", "--has_vr", help = "Does the table have VR? The VR column is assumed to be titled HRV", action = "store_true")
parser.add_argument("-m", "--has_met", help = "Does the table have metallicity/alpha? The columns are assumed to be titled Alpha_c and Met_N_K", action = "store_true")
parser.add_argument("-t", "--threshold", help = "-1 * Log-likelihood threshold for component membership", default = 9.0, type = float)
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

t = gfc.io.read_csv(args.data_file)

if args.verbose:
    print "Finished loading data"

amps_xd = np.load(args.xd_results_folder+"amplitudes.npy")
means_xd = np.load(args.xd_results_folder+"means.npy")
covs_xd = np.load(args.xd_results_folder+"covariances.npy")

if args.verbose:
    print "Finished loading XD parameters"

gfc.add_rad(t)
if args.has_vr:
    gfc.matrix.add_w(t, v_r_col = "HRV")
else:
    gfc.matrix.add_w(t)
gfc.matrix.add_A(t)
gfc.matrix.add_R(t)
gfc.matrix.add_UVW(t)

if args.verbose:
    print "Calculated U, V, W"

Ls = gfc.pdf.loglikelihood_many_multiPDF(means_xd, covs_xd, t["UVW_vec"])

if args.verbose:
    print "Calculated likelihoods"

met, alpha = (t["Met_N_K"], t["Alpha_c"]) if args.has_met else (None, None)

gfc.gplot.moving_group(t["ra"], t["dec"], t["l"], t["b"], t["U"], t["V"], t["W"], range(len(t)), [], met = met, alpha = alpha, bins = 100, saveto = "{f}/all_mg.png".format(f = args.save_folder)) 
for n, (a, m, c) in enumerate(zip(amps_xd, means_xd, covs_xd)):
    indices_in = np.where(Ls[:,n] > -args.threshold)[0]
    indices_ex = np.where(Ls[:,n] <= -args.threshold)[0]
    if args.verbose:
        print "{n:02d}: {l}".format(n=n+1, l=len(indices_in))
    gfc.gplot.moving_group(t["ra"], t["dec"], t["l"], t["b"], t["U"], t["V"], t["W"], indices_in, indices_ex, met = met, alpha = alpha, a = a, m = m, c = c, bins = 100, saveto = "{f}/{n:02d}_mg.png".format(f = args.save_folder, n = n+1)) 

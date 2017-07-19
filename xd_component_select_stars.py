import gfc
import numpy as np

from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("data_file", help = "File that contains the stellar data")
parser.add_argument("xd_results_folder", help = "Folder that contains the XD results")
parser.add_argument("save_folder", help = "Folder in which to save results")
parser.add_argument("-v", "--verbose", action = "store_true")
parser.add_argument("-t", "--threshold", help = "Log-likelihood threshold for component membership", default = -9)
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

#PDFs = [gfc.pdf.multivariate(m, c, a) for a, m, c in zip(amps_xd, means_xd, covs_xd)]
#Ls = gfc.pdf.likelihood_many(PDFs, t["UVW_vec"])

Ls = gfc.pdf.loglikelihood_many_multiPDF(means_xd, covs_xd, t["UVW_vec"])


if args.verbose:
    print "Calculated likelihoods"

print Ls.shape

print Ls.max(axis = 0)
print Ls.min(axis = 0)
print np.mean(Ls, axis = 0)
print np.median(Ls, axis = 0)

t_subs = [t[Ls[:,n] > -9] for n in range(len(means_xd))]
for n, t_ in enumerate(t_subs):
    print n, len(t_)
    gfc.gplot.moving_group(t_["ra"], t_["dec"], t_["l"], t_["b"], t_["U"], t_["V"], t_["W"], amps_xd[n], means_xd[n], covs_xd[n], bins = 100, saveto = "{0}_mg.png".format(n+1)) 

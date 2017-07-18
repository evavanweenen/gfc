import gfc
import numpy as np

from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("data_file", help = "File that contains the stellar data")
parser.add_argument("xd_results_folder", help = "Folder that contains the XD results")
parser.add_argument("save_folder", help = "Folder in which to save results")
parser.add_argument("-v", "--verbose", action = "store_true")
parser.add_argument("-n", help = "Which component to plot", type = int, default = 0)
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

PDFs = [gfc.pdf.multivariate(m, c, a) for a, m, c in zip(amps_xd, means_xd, covs_xd)]
Ls = gfc.pdf.likelihood_many(PDFs, t["UVW_vec"])
    
if args.verbose:
    print "Calculated likelihoods"

for n in range(19):
    t.add_column(gfc.table.Column(data = Ls[:, n], name = "L"))
    t.sort("L")
    t.reverse()
    gfc.gplot.plt.hist2d(t["U"][:5000], t["W"][:5000], range = ((-120, 120), (-70, 70)), bins = 100)
    gfc.gplot.plt.xlim(-120,120)
    gfc.gplot.plt.ylim(-70,70)
    gfc.gplot.draw_PDF_ellipse(gfc.gplot.plt.gca(), amps_xd[n], means_xd[n], covs_xd[n], "xz")
    gfc.gplot.plt.savefig(args.save_folder + "stars_{n}.png".format(n = n))
    gfc.gplot.plt.close()
    t.remove_column("L")
#for N in range(19):
#    Ls[:,N] = np.log(Ls[:,N])
#    print "{N}: min {min:.1f} ; median {median:.1f} ; 95% {ten:.1f} ; max {max:.1f} ".format(N = N, max = np.nanmax(Ls[:,N]), min = np.nanmin(Ls[:,N]), median = np.median(Ls[:,N]), ten = np.percentile(Ls[:,N], 95))

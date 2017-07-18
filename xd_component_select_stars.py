import gfc
import numpy as np

from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("data_file", help = "File that contains the stellar data")
parser.add_argument("xd_results_folder", help = "Folder that contains the XD results")
parser.add_argument("save_folder", help = "Folder in which to save results")
parser.add_argument("-v", "--verbose", action = "store_true")
parser.add_argument("-n", help = "Which component to plot", type = int, default = 0)
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

t.add_column(gfc.table.Column(data = Ls[:, args.n], name = "L"))
t.sort("L")
t.reverse()
gfc.gplot.plt.scatter(t["U"][:5000], t["W"][:5000])
gfc.gplot.plt.xlim(-120,120)
gfc.gplot.plt.ylim(-70,70)
gfc.gplot.draw_PDF_ellipse(gfc.gplot.plt.gca(), amps_xd[args.n], means_xd[args.n], covs_xd[args.n], "xz")
gfc.gplot.plt.savefig(args.save_folder + "stars.png")

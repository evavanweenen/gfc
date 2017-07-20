import matplotlib
matplotlib.use('Agg') 
import gfc
import numpy as np
from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("data_file", help = "File that contains the stellar data")
parser.add_argument("xd_results_folder", help = "Folder that contains the XD results")
parser.add_argument("save_folder", help = "Folder in which to save results")
parser.add_argument("-t", "--threshold", type = float, help = "Threshold angle to select stars", default = 2.5)
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

time_before_reading = gfc.time()
t = gfc.io.read_csv(args.data_file)
time_after_reading = gfc.time()
if args.verbose:
    print "Finished loading data in {0:.1f} seconds".format(time_after_reading - time_before_reading)

gfc.add_rad(t)

amps, means, covs = gfc.io.load_PDFs(args.xd_results_folder)
if args.verbose:
    print "Finished loading XD parameters"

time_before_dot = gfc.time()
mean_coords = gfc.matrix.mean_to_coords_many(means)
dot_products = gfc.matrix.mean_coords_dot_stars(mean_coords, t["ra_rad"], t["dec_rad"])
time_after_dot = gfc.time()
if args.verbose:
    print "Calculated dot products in {0:.1f} seconds".format(time_after_dot - time_before_dot)

crit = 1. - np.cos(args.threshold * np.pi / 180.)

selections = [np.where(np.abs(d) < crit)[0] for d in dot_products]

for i, (m, s) in enumerate(zip(mean_coords, selections)):
    pml, pmb = gfc.matrix.ICRS_to_galactic.transformProperMotions(t["ra_rad"][s], t["dec_rad"][s], t["pmra"][s], t["pmdec"][s])
    fig, axs = gfc.gplot.plt.subplots(2, tight_layout = True)
    axs[0].scatter(m[0]*180/np.pi, m[1]*180/np.pi, s = 20, c='k')
    #gfc.gplot.plt.scatter(t["ra_rad"][s], t["dec_rad"][s], s = 1, color='r')
    gfc.gplot.density_ax(axs[0], t["ra"][s], t["dec"][s], r = ((0, 360), (-90, 90)), bins = (75, 40))
    axs[1].hist(pmb, bins = 25, range = (-150, 150))
    axs[1].set_xlim(-120, 120)
    axs[1].set_xlabel("$\mu(b)$")
    fig.savefig("{f}/ring_{n:02d}.png".format(f = args.save_folder, n = i))
    gfc.gplot.plt.close()
    print "{n:02d}: {N:04d} stars".format(n = i + 1, N = len(s))

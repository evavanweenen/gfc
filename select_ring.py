import gfc
import numpy as np

from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("data_file", help = "File that contains the stellar data")
parser.add_argument("xd_results_folder", help = "Folder that contains the XD results")
parser.add_argument("save_folder", help = "Folder in which to save results")
parser.add_argument("-t", "--threshold", type = float, help = "Threshold angle to select stars", default = 1.)
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
    #gfc.gplot.plt.scatter(m[0], m[1], s = 20, c='k')
    #gfc.gplot.plt.scatter(t["ra_rad"][s], t["dec_rad"][s], s = 1, color='r')
    gfc.gplot.density(t["ra"][s], t["dec"][s], r = ((0, 360), (-90, 90)), bins = 100, saveto = "{f}/ring_{n}.png".format(f = args.save_folder, n = i))
    print "{n:02d}: {N:04d} stars".format(n = i + 1, N = len(s))

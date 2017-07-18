"""
Make HR diagrams of the data with and without the main sequence
"""

import gfc
import matplotlib as mpl

from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("data_file", help = "File with data to be plotted")
parser.add_argument("-s", "--saveto", help = "Location to save figure to", default = "CMD.png")
parser.add_argument("-x", "--x_col", help = "Column to plot on x-axis", default = "g_min_ks")
parser.add_argument("-y", "--y_col", help = "Column to plot on y-axis", default = "g_mag_abs")
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

t = gfc.io.read_csv(args.data_file)
if args.verbose:
    print "Finished reading" 
gfc.gplot.density(t[args.x_col], t[args.y_col], saveto = args.saveto, xlabel = "$G - K$", ylabel = "$G$", flip = "y", xlim = (-0.5, 4.0), ylim = (11, -1.8), cb=True, bins = 175, xticks = range(0, 5), title = "Final TGAS sample", histkwargs = {"vmin":1, "vmax":1200})
if args.verbose:
    print "Finished plotting"

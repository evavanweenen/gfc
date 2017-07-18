import gfc
import numpy as np

from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("data_file", help = "File that contains the stellar data")
parser.add_argument("save_folder", help = "Folder in which to save results")
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

t = gfc.io.read_csv(args.data_file)

if args.verbose:
    print "Finished loading data"

gfc.tgas.add_rad(t)
gfc.tgas.add_w(t)
gfc.tgas.add_A(t)
gfc.tgas.add_R(t)
gfc.tgas.add_UVW(t)

if args.verbose:
    print "Calculated U, V, W"

gfc.gplot.density(t["U"], t["W"], bins = 200, r = ((-120, 120), (-70, 70)))

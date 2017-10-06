import numpy as np
import gfc
from astropy.table import Table
from gfc import ArgumentParser

parser = ArgumentParser()
parser.add_argument("save_folder", help = "Folder in which results will be saved")
parser.add_argument("--L", help = "How many stars to simulate", default=1e6, type = int)
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

if args.verbose:
    print "Finished parsing arguments"


"""
Simulate data
"""
args.L = int(args.L)
arr=np.empty([args.L, 10], dtype=float)
t=Table(arr, names=('ra', 'dec', 'ra_error', 'dec_error', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'pmra_error', 'pmdec_error'))
t['ra'][:] = np.random.normal(180, 0.341*180, args.L)
t['dec'][:] = np.random.normal(180, 0.341*180, args.L)
t['ra_error'][:] = np.random.normal(.1, .05, args.L)
t['dec_error'][:] = np.random.normal(.1, .03, args.L)
t['parallax'][:] = np.random.normal(5, .1, args.L)
t['parallax_error'][:] = np.random.normal(.1, .05, args.L)
t['pmra'][:] = np.random.normal(30, 10, args.L)
t['pmdec'][:] = np.random.normal(30, 10, args.L)
t['pmra_error'][:] = np.random.normal(.1, .05, args.L)
t['pmdec_error'][:] = np.random.normal(.1, .03, args.L)

print "t", t

gfc.io.write_csv(t, save_folder)

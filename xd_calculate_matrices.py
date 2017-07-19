"""
Calculate several matrices used in XD.
Run this script to do it once so it doesn't have to be re-done every time (can be quite computationally heavy).
"""

import gfc
from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("data_file", help = "File containing TGAS data")
parser.add_argument("save_folder", help = "Folder in which results will be saved")
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

if args.verbose:
    print "Now reading data..."
time_before_reading = gfc.time()
t = gfc.io.read_csv(args.data_file) # ASSUME CSV -- todo: check
if args.verbose:
    print "{0} lines in data".format(len(t))
    print "Data contains the following keys:"
    print t.keys()

gfc.tgas.add_rad(t)
gfc.remove_unused_columns(t)
time_after_reading = gfc.time()
if args.verbose:
    print "Read TGAS data in {0:.1f} seconds".format(time_after_reading - time_before_reading)

time_before_matrices = gfc.time()
gfc.tgas.add_w(t)
gfc.tgas.add_A(t)
gfc.tgas.add_R(t)
gfc.tgas.add_C(t)
gfc.tgas.add_Q(t)
gfc.tgas.add_S(t)
time_after_matrices = gfc.time()
if args.verbose:
    print "Added w, A, R, C, Q, S in {0:.1f} seconds".format(time_after_matrices - time_before_matrices)

time_before_writing = gfc.time()
gfc.io.write_table_with_separate_arrays(t, saveto_folder = args.save_folder, verbose = args.verbose)
time_after_writing = gfc.time()
if args.verbose:
    print "Written in {0:.1f} seconds".format(time_after_writing - time_before_writing)

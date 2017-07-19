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
if args.verbose:
    print "Added w"

gfc.tgas.add_A(t)
if args.verbose:
    print "Added A"

gfc.tgas.add_R(t)
if args.verbose:
    print "Added R and R^-1"

gfc.tgas.add_C(t)
if args.verbose:
    print "Added C"

gfc.tgas.add_Q(t)
if args.verbose:
    print "Added Q"

gfc.tgas.add_S(t)
if args.verbose:
    print "Added S"
time_after_matrices = gfc.time()
if args.verbose:
    print "Calculated matrices in {0:.1f} seconds".format(time_after_matrices - time_before_matrices)

time_before_writing = gfc.time()
gfc.io.write_table_with_separate_arrays(t, saveto_folder = args.save_folder, verbose = args.verbose)
time_after_writing = gfc.time()
if args.verbose:
    print "Written in {0:.1f} seconds".format(time_after_writing - time_before_writing)

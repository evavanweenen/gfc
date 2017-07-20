import gfc

from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("tmass_file", help = "File containing the TGAS/2MASS table")
parser.add_argument("rave_file", help = "File containing the RAVE table")
parser.add_argument("save_to", help = "Location to save cross-matched table to")
parser.add_argument("-v", "--verbose", action = "store_true")
parser.add_argument("-r", "--remove_columns", action = "store_true")
parser.add_argument("-a", "--add_velocities", action = "store_true")
args = parser.parse_args()

if args.verbose:
    print "Now loading TGAS/2MASS table"
t = gfc.io.read_csv(args.tmass_file)
t.rename_column("tycho2_id", "ID_TYCHO2")
t.rename_column("source_id", "ID")
original_length = len(t)
t.remove_rows(gfc.np.where(t["ID_TYCHO2"].mask)[0])
new_length = len(t)
if args.verbose:
    print "Loaded TGAS/2MASS table; removed {0}/{1} rows without a TYCHO2 ID".format(original_length - new_length, original_length)
if args.verbose:
    print "Now loading RAVE table"
rave = gfc.io.read_csv(args.rave_file)
original_length = len(rave)
rave.remove_rows(gfc.np.where(rave["ID_TYCHO2"].mask)[0])
new_length = len(rave)
if args.verbose:
    print "Loaded RAVE table; removed {0}/{1} rows without a TYCHO2 ID".format(original_length - new_length, original_length)
if args.verbose:
    print "Now performing cross-match"
j = gfc.table.join(rave, t, keys="ID_TYCHO2")
del rave, t
if args.verbose:
    print "Cross-match done: {0} rows included".format(len(j))
if args.remove_columns:
    gfc.remove_unused_columns(j)
if args.add_velocities:
    gfc.add_rad(j)
    gfc.matrix.add_w(j, v_r_col = "HRV", components = False)
    gfc.matrix.add_A(j)
    gfc.matrix.add_R(j)
    gfc.matrix.add_UVW(j, vector = False)
    j.remove_columns(("R", "A", "R^-1", "w_vec", "ra_rad", "dec_rad", "ra_rad_error", "dec_rad_error"))
gfc.io.write_csv(j, args.save_to)
if args.verbose:
    print "Cross-matched table written"

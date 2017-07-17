import gfc
from sys import argv

TMASS = argv[1] # folder?
RAVE  = argv[2] # CSV

print "Now loading 2MASS table"
t = gfc.io.load_table_with_separate_arrays(saveto_folder="tgas/2MASS/results/")
t.remove_columns(("R", "Q", "S", "A", "R^-1", "w1", "w2", "w3"))
t.rename_column("tycho2_id", "ID_TYCHO2")
t.rename_column("source_id", "ID")
t.remove_rows(gfc.np.where(t["ID_TYCHO2"].mask)[0])
print "Table loaded"
print "Now loading RAVE table"
rave = gfc.io.read_csv(RAVE)
rave.remove_rows(gfc.np.where(rave["ID_TYCHO2"].mask)[0])
print "Table loaded"

joinedtable = gfc.table.join(rave, t, keys="ID_TYCHO2")

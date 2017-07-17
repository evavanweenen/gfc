"""
Plot a histogram of local velocity space using TGAS and RAVE data
"""

import gaia_fc as g
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams.update({"font.size": 25, "figure.figsize": (10, 10)})
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec


print "imported gaia_fc"

try:
    assert t
except:
    print "Now loading 2MASS table"
    t = g.io.load_table_with_separate_arrays(saveto_folder="tgas/2MASS/results/")
    t.remove_columns(("R", "Q", "S", "A", "R^-1", "w1", "w2", "w3"))
    t.rename_column("tycho2_id", "ID_TYCHO2")
    t.rename_column("source_id", "ID")
    t.remove_rows(g.np.where(t["ID_TYCHO2"].mask)[0])
    print "Table loaded"
    print "Now loading RAVE table"
    rave = g.io.read_csv("tgas/RAVE/RAVE_TGAS_1000.csv")
    rave.remove_rows(g.np.where(rave["ID_TYCHO2"].mask)[0])
    print "Table loaded"

    joinedtable = g.table.join(rave, t, keys="ID_TYCHO2")
    del t
    del rave
    t = joinedtable
try:
    assert summary
except:
    summary = g.io.read("tgas/summary_rave.txt", format="fixed_width")
    t = g.table.join(t, summary, keys="ID")

if "U" not in t.keys():
    t["Glon"] = g.radians(t["Glon"])
    t["Glat"] = g.radians(t["Glat"])
    t["RA_TGAS"] = g.radians(t["RA_TGAS"])
    t["DE_TGAS"] = g.radians(t["DE_TGAS"])

    t["pmGlon"], t["pmGlat"] = zip(*[g.ICRS_to_galactic.transformProperMotions(p,th,mp,mth) for p,th,mp,mth in zip(t["RA_TGAS"], t["DE_TGAS"], t["pmRA_TGAS"], t["pmDE_TGAS"])])

    v_totR = g.total_velocity(t["pmRA_TGAS"], t["pmDE_TGAS"], t["parallax_TGAS"], t["HRV_1"])
    v_totR.name = "v_totR"
    t.add_column(v_totR)

    x,y,z,UR,VR,WR = g.tophase(t["Glon"], t["Glat"], t["parallax_TGAS"], t["pmGlon"], t["pmGlat"], t["HRV_1"])
    x.name, y.name, z.name, UR.name, VR.name, WR.name = "x", "y", "z", "U", "V", "W"
    t.add_columns((x,y,z,UR,VR,WR))
    print "calculated observed velocities"

g.gplot.density(t["U"], t["V"], xlabel="$U$ (km/s)", ylabel="$V$ (km/s)", saveto="UVplane.png", r=((-130,120),(-120,60)), bins = 125, yticks=(-100,-50,0,50))

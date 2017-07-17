"""
Calculate S^2, U, V, W velocities as functions of B-V colours to:
- look for Parenago's discontinuity
- determine the motion of the Sun
"""

import gaia_fc as g

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

if "v_totM" not in t.keys():
    t["Glon"] = g.radians(t["Glon"])
    t["Glat"] = g.radians(t["Glat"])
    t["RA_TGAS"] = g.radians(t["RA_TGAS"])
    t["DE_TGAS"] = g.radians(t["DE_TGAS"])

    t["pmGlon"], t["pmGlat"] = zip(*[g.ICRS_to_galactic.transformProperMotions(p,th,mp,mth) for p,th,mp,mth in zip(t["RA_TGAS"], t["DE_TGAS"], t["pmRA_TGAS"], t["pmDE_TGAS"])])

    v_totR = g.total_velocity(t["pmRA_TGAS"], t["pmDE_TGAS"], t["parallax_TGAS"], t["HRV_1"])
    v_totR.name = "v_totR"
    t.add_column(v_totR)

    x,y,z,UR,VR,WR = g.tophase(t["Glon"], t["Glat"], t["parallax_TGAS"], t["pmGlon"], t["pmGlat"], t["HRV_1"])
    x.name, y.name, z.name, UR.name, VR.name, WR.name = "x", "y", "z", "UR", "VR", "WR"
    t.add_columns((x,y,z,UR,VR,WR))
    print "calculated observed velocities"
    
    v_totM = g.total_velocity(t["pmRA_TGAS"], t["pmDE_TGAS"], t["parallax_TGAS"], t["best_c"])
    v_totM.name = "v_totM"
    t.add_column(v_totM)

    x,y,z,UM,VM,WM = g.tophase(t["Glon"], t["Glat"], t["parallax_TGAS"], t["pmGlon"], t["pmGlat"], t["best_c"])
    x.name, y.name, z.name, UM.name, VM.name, WM.name = "x", "y", "z", "UM", "VM", "WM"
    t.add_columns((UM,VM,WM))
    print "calculated model velocities"

B = (1000,0.04); range_binover=(-0.2, 1.5)
(meansUR, sigmasUR), binsUR, highestUR, sigmaBV = g.mapping.map_over_bins((g.np.mean, g.mean_sigma), t["BminV"], B, -t["UR"], range_binover=range_binover)
(meansVR, sigmasVR), binsVR, highestVR, sigmaBV = g.mapping.map_over_bins((g.np.mean, g.mean_sigma), t["BminV"], B, -t["VR"], range_binover=range_binover)
(meansWR, sigmasWR), binsWR, highestWR, sigmaBV = g.mapping.map_over_bins((g.np.mean, g.mean_sigma), t["BminV"], B, -t["WR"], range_binover=range_binover)
(S2R, sigmaS2R), binsS2R, highestS2R, sigmaBV = g.mapping.map_over_bins((g.S2, g.sigma_S2), t["BminV"], B, t["RA_TGAS", "DE_TGAS", "UR", "VR", "WR"], range_binover=range_binover)
UR = [binsUR, meansUR, sigmasUR, highestUR]
VR = [binsVR, meansVR, sigmasVR, highestVR]
WR = [binsWR, meansWR, sigmasWR, highestWR]
S2R_ = [binsS2R, g.np.sqrt(S2R), None, highestS2R]

(meansUM, sigmasUM), binsUM, highestUM, sigmaBV = g.mapping.map_over_bins((g.np.mean, g.mean_sigma), t["BminV"], B, -t["UM"], range_binover=range_binover)
(meansVM, sigmasVM), binsVM, highestVM, sigmaBV = g.mapping.map_over_bins((g.np.mean, g.mean_sigma), t["BminV"], B, -t["VM"], range_binover=range_binover)
(meansWM, sigmasWM), binsWM, highestWM, sigmaBV = g.mapping.map_over_bins((g.np.mean, g.mean_sigma), t["BminV"], B, -t["WM"], range_binover=range_binover)
(S2M, sigmaS2M), binsS2M, highestS2M, sigmaBV = g.mapping.map_over_bins((g.S2, g.sigma_S2), t["BminV"], B, t["RA_TGAS", "DE_TGAS", "UM", "VM", "WM"], range_binover=range_binover)
UM = [binsUM, meansUM, sigmasUM, highestUM]
VM = [binsVM, meansVM, sigmasVM, highestVM]
WM = [binsWM, meansWM, sigmasWM, highestWM]
S2M_ = [binsS2M, g.np.sqrt(S2M), None, highestS2M]
g.gplot.BV_SUVW(UR, VR, WR, S2R_, UM, VM, WM, S2M_, saveto="BV-UVWS.png")
g.gplot.BV_SV(VR, S2R_, VM, S2M_, saveto="BV-VS.png")

tblue = t[t["BminV"] > 0]
UmeanR = -tblue["UR"].mean() ; UmeanerrR = g.mean_sigma(tblue["UR"])
WmeanR = -tblue["WR"].mean() ; WmeanerrR = g.mean_sigma(tblue["WR"])
VcoeR, VcovR = g.np.polyfit(S2R[binsS2R > 0], meansVR[binsS2R > 0], 1, cov=True)
Vat0R = VcoeR[1] ; Vat0errR = g.np.sqrt(VcovR[1,1])
print "Solar velocity (RAVE):"
print "U: {0:.2f} +- {1:.2f}".format(UmeanR, UmeanerrR)
print "V: {0:.2f} +- {1:.2f}".format(Vat0R, Vat0errR)
print "W: {0:.2f} +- {1:.2f}".format(WmeanR, WmeanerrR)
UmeanM = -tblue["UM"].mean() ; UmeanerrM = g.mean_sigma(tblue["UM"])
WmeanM = -tblue["WM"].mean() ; WmeanerrM = g.mean_sigma(tblue["WM"])
VcoeM, VcovM = g.np.polyfit(S2M[binsS2M > 0], meansVM[binsS2M > 0], 1, cov=True)
Vat0M = VcoeM[1] ; Vat0errM = g.np.sqrt(VcovM[1,1])
print "Solar velocity (XD):"
print "U: {0:.2f} +- {1:.2f}".format(UmeanM, UmeanerrM)
print "V: {0:.2f} +- {1:.2f}".format(Vat0M, Vat0errM)
print "W: {0:.2f} +- {1:.2f}".format(WmeanM, WmeanerrM)
g.gplot.S2UVW(UR, VR, WR, S2R_, UM, VM, WM, S2M_, UmeanR, VcoeR, WmeanR, UmeanM, VcoeM, WmeanM, saveto="S2-UVW.png")
g.gplot.S2V(VR, S2R_, VM, S2M_, VcoeR, VcoeM, saveto="S2-V.png")

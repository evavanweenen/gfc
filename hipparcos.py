import numpy as np
from gaia_fc import general as gen
from astropy.io import votable as vot, ascii
from astropy import table
from matplotlib import pyplot as plt

ex = 131

p = {'c': 'k', 's': 1}
extent = (-ex, ex, -ex, ex)

level_min = 0.0000001
level_max = 0.01
levels = np.logspace(np.log10(level_min), np.log10(level_max), 15)

def read_data_cols(t):
    x = [("plx", "mas", r"$\bar\omega$"),
         ("ra_rad", "rad", r"$\alpha$"),
         ("de_rad", "rad", r"$\delta$"),
         ("pm_ra", "mas/yr", r"$\mu_\alpha$"),
         ("pm_de", "mas/yr", r"$\mu_\delta$"),
         ("hp_mag", "mag", r"$H_p$"),
         ("m", "mag", r"$H_p$"),
         ("b_v", "mag", r"$B - V$")]

    for col, unit, label in x:
        t[col].unit = unit
        t[col].axis_label = label

def read_data_as_table(loc = "hip/hipparcos.vot"):
    data = vot.parse(loc, unit_format="vounit")
    t = data.get_first_table().to_table()
    del data # clear memory
    read_data_cols(t)
    return t

def read_multiple_data_files_as_tables(*locs):
    tables = [read_data_as_table(loc) for loc in locs]
    t = table.vstack(tables)
    t = table.unique(t, "hip")
    tycho = ascii.read("hip/hip_after_tycho.txt", format="fixed_width")
    t = t[[j for j,h in enumerate(t["hip"]) if h in tycho["hip"]]]
    read_data_cols(t)
    t.remove_column("hipparcos_bovy2_oid")
    return t

def read():
    table_south = read_data_as_table("hip/hip_south.vot")
    table_mag = ascii.read("hip/hip_mag.txt", format="fixed_width")
    t = table.vstack((table_mag, table_south))
    t = table.unique(t, "hip")
    read_data_cols(t)
    return t

def plot_columns(t, saveto, *cols, **kw):
    cols_not_in_t = [col for col in cols if col not in t.keys()]
    if cols_not_in_t:
        raise gen.GFC_Exception("\n\ngaia_fc.hipparcos.plot_columns: The following columns are missing in the table: {0}".format(cols_not_in_t))
    for j, p1 in enumerate(cols):
        for p2 in cols[j+1:]:
            plt.scatter(t[p1], t[p2], **p)
            plt.xlabel(t[p1].axis_label + " (" + t[p1].unit.to_string() + ")")
            plt.ylabel(t[p2].axis_label + " (" + t[p2].unit.to_string() + ")")
            plt.title("Hipparcos: {0} vs {1}".format(p1, p2))
            plt.axis("tight")
            if kw:
                plt.setp(plt.gca(), **kw)
            plt.tight_layout()
            plt.savefig(saveto.format(p1, p2))
            plt.close()

def plot_pdf_and_v_data(t, saveto, PDF, title="{0} - {1}"):
    "hard coded to work with v_x, v_y, v_z"

    cols_not_in_t = [col for col in ("v_x", "v_y", "v_z") if col not in t.keys()]
    if cols_not_in_t:
        raise gen.GFC_Exception("\n\ngaia_fc.hipparcos.plot_pdf_and_v_data: The following columns are missing in the table: {0}".format(cols_not_in_t))
    if PDF.ndim != 3:
        raise gen.GFC_Exception("\n\ngaia_fc.hipparcos.plot_pdf_and_v_data: The provided PDF has dimensions {0} -- expected 3 dimensions".format(PDF.ndim))

    for i, pair, xl, yl in zip([0, 1, 2], [("v_y", "v_z"), ("v_x", "v_z"), ("v_x", "v_y")], [(-130, 120), (-130, 120), (-130, 120)], [(-70, 70), (-70, 70), (-120, 60)]):
        p1, p2 = pair
        pair_int = [0, 1, 2]
        pair_int.remove(i)
        pdfsumi = PDF.sum(axis=i).T
        plt.scatter(t[p1], t[p2], s = 1, c = "0.2", edgecolors="face")
        plt.contour(pdfsumi, extent=extent, linewidths=3, levels=levels, colors='k')
        plt.xlabel(t[p1].axis_label + " (" + t[p1].unit.to_string() + ")")
        plt.ylabel(t[p2].axis_label + " (" + t[p2].unit.to_string() + ")")
        plt.title(title.format(p1, p2))
        plt.xlim(xl)
        plt.ylim(yl)
        plt.tight_layout()
        plt.savefig(saveto.format(p1, p2))
        plt.close()

def add_lb(t):
    tran = gen.coords.CoordinateTransformation(gen.coords.Transformations.ICRS2GAL)
    l, b = tran.transformSkyCoordinates(t["ra_rad"], t["de_rad"])
    t.add_column(table.Column(data = l, name = "l", unit = "rad"))
    t.add_column(table.Column(data = b, name = "b", unit = "rad"))

def add_w(t):
    ws = gen.map_np(gen.w, t["de_rad"], t["plx"], t["pm_ra"], t["pm_de"], np.zeros(len(t)))
    t.add_column(table.Column(data = ws[:, 0, 0], name = "w1", unit = "km / s"))
    t.add_column(table.Column(data = ws[:, 1, 0], name = "w2", unit = "1 / yr"))
    t.add_column(table.Column(data = ws[:, 2, 0], name = "w3", unit = "1 / yr"))
    t["w1"].axis_label = "$w_1$"
    t["w2"].axis_label = "$w_2$"
    t["w3"].axis_label = "$w_3$"

def add_A(t):
    As = gen.map_np(gen.A, t["ra_rad"], t["de_rad"])
    t.add_column(table.Column(data = As, name = "A"))

def add_R(t):
    R_invs = gen.pmap_np(gen.R_inv, t["A"])
    t.add_column(table.Column(data = R_invs, name = "R^-1"))
    Rs = np.linalg.inv(R_invs)
    t.add_column(table.Column(data = Rs, name = "R"))

def C(U):
    uinv = np.linalg.inv(U)
    return uinv.dot(uinv.T)
def C_row_old(row):
    U1 = np.array([row["u1"], row["u2"], row["u4"], row["u7"], row["u11"]])
    U2 = np.array([0., row["u3"], row["u5"], row["u8"], row["u12"]])
    U3 = np.array([0., 0., row["u6"], row["u9"], row["u13"]])
    U4 = np.array([0., 0., 0., row["u10"], row["u14"]])
    U5 = np.array([0., 0., 0., 0., row["u15"]])
    U = np.vstack((U1, U2, U3, U4, U5))
    C_here = C(U)

    return C_here
def C_row(row):
    nu = row["ntr"] - 5
    sigma = np.array([row["e_ra_rad"], row["e_de_rad"], row["e_plx"], row["e_pm_ra"], row["e_pm_de"]])
    nu29 = 2./(9. + nu)
    Q = nu * np.power((np.sqrt(nu29) * row["f2"] + 1. - nu29), 3)
    u = np.sqrt(Q / nu)
    U = np.array([[row["u1"], row["u2"], row["u4"], row["u7"], row["u11"]],
                  [0.,        row["u3"], row["u5"], row["u8"], row["u12"]],
                  [0.,        0.,        row["u6"], row["u9"], row["u13"]],
                  [0.,        0.,        0.,        row["u10"],row["u14"]],
                  [0.,        0.,        0.,        0.,        row["u15"]]])
    for i, s in enumerate(sigma):
        U[i, i] *= u/s

    U_prod = U.T.dot(U)
    U_prod[[0, 1], :] *= gen.radiantomas # RA and Dec rows    -- colon for clarity
    U_prod[:, [0, 1]] *= gen.radiantomas # RA and Dec columns -- double for RARA/RADec/DecDec terms is intentional

    C = np.linalg.inv(U_prod) * u**2. # n.b. second u is the scalar

    return C
def add_C(t):
    Cs = gen.map_np(C_row, t) # pmap crashes
    t.add_column(table.Column(data = Cs, name = "C"))

def Q_row(row):
    return gen.Q(row["ra_rad"], row["de_rad"], row["plx"], row["pm_ra"], row["pm_de"])
def add_Q(t):
    Qs = gen.map_np(gen.Q_star, t["plx"], t["pm_ra"], t["pm_de"]) # no pmap
    t.add_column(table.Column(data = Qs, name = "Q"))

def add_S(t):
    Ss = gen.map_np(gen.S, t["C"], t["Q"])
    t.add_column(table.Column(data = Ss, name = "S"))

def add_v_proj(t):
    ws = gen.XD_arr(t, "w1", "w2", "w3")
    vs = gen.map_np(gen.v_proj, ws, t["R^-1"])

    t.add_column(table.Column(data = vs[:, 0], name = "v_x", unit = "km / s")) ; t["v_x"].axis_label = "$v_x$"
    t.add_column(table.Column(data = vs[:, 1], name = "v_y", unit = "km / s")) ; t["v_y"].axis_label = "$v_y$"
    t.add_column(table.Column(data = vs[:, 2], name = "v_z", unit = "km / s")) ; t["v_z"].axis_label = "$v_z$"

def rhat_row(row):
    return gen.rhat(row["ra_rad"], row["de_rad"])
def add_rhat(t):
    rhats = gen.map_np(rhat_row, t)
    t.add_column(table.Column(data = rhats, name = "r_hat"))

def add_Rr_and_Rt(t):
    Rrs = gen.map_np(gen.R_r, t["r_hat"])
    Rts = np.identity(3) - Rrs

    t.add_column(table.Column(data = Rrs, name = "R_r"))
    t.add_column(table.Column(data = Rts, name = "R_t"))

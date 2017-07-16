import numpy as np
from astropy import table
from matplotlib import pyplot as plt
from gaia_fc import general as gen

p = {"alpha": 0.5, "cmap": plt.cm.rainbow}

extent = (-200, 200, -200, 200)

level_min = 0.0000001
level_max = 0.01
levels = np.logspace(np.log10(level_min), np.log10(level_max), 15)

def generate_velocities(mean_, cov_, nr):
    mean = np.array(mean_)
    cov = gen.covariant_array(cov_)
    velocities = np.random.multivariate_normal(mean, cov, size=nr)
    err = np.ones_like(velocities)
    velocities = np.hstack((velocities, err))
    return velocities

def generate_positions(N = 1e5, R = 100.):
    x = np.random.uniform(-R, R, size=2*N)
    y = np.random.uniform(-R, R, size=2*N)
    z = np.random.uniform(-R, R, size=2*N)

    dists = np.sqrt(x**2. + y**2. + z**2.)

    x = x[dists <= R]
    y = y[dists <= R]
    z = z[dists <= R]
    a = np.column_stack((x, y, z))
    np.random.shuffle(a)
    a = a[:N]
    a = np.hstack((a, np.ones_like(a)))

    return a

def generate_star_group(N, R, v_mean, v_cov, label):
    vel = generate_velocities(v_mean, v_cov, N)
    pos = generate_positions(N, R)
    lab = np.tile(label, N)

    stars = np.hstack((pos, vel, lab[:, np.newaxis]))
    return stars

def star_group_to_table(group):
    stars = table.Table(data = group, names=["x", "y", "z", "e_x", "e_y", "e_z", "v_x", "v_y", "v_z", "e_v_x", "e_v_y", "e_v_z", "group"], dtype = ["f8"]*12+["i1"])
    for ind in ("x", "y", "z", "e_x", "e_y", "e_z"):
        stars[ind].unit = "pc"
    for ind in ("v_x", "v_y", "v_z", "e_v_x", "e_v_y", "e_v_z"):
        stars[ind].unit = "km/s"
    for ind in ("x", "y", "z", "e_x", "e_y", "e_z", "v_x", "v_y", "v_z"):
        stars[ind].axis_label = "${0}$".format(ind)
    for ind in ("e_v_x", "e_v_y", "e_v_z"):
        stars[ind].axis_label = "$" + r"e_{" + ind[2:] + r"}$"
    return stars

def generate_stars(*groups):
    """
    group:
        [N, R, v_tuple, cov_tuple, label]
    """
    star_groups = [generate_star_group(*group) for group in groups]
    all_stars = np.concatenate(star_groups)
    np.random.shuffle(all_stars)

    t = star_group_to_table(all_stars)

    return t

def add_astrometry(t):
    cols_already_in_t = [col for col in ("phi", "theta", "plx", "muphi*", "mutheta", "v_r") if col in t.keys()]
    if cols_already_in_t:
        raise gen.GFC_Exception("\n\ngaia_fc.simulated.add_astrometry: The following columns are already in the table: {0}".format(cols_already_in_t))

    as_astrometry = np.array(gen.toastro(t["x"], t["y"], t["z"], t["v_x"], t["v_y"], t["v_z"]))
    t.add_column(table.Column(data = as_astrometry[0], name="phi", unit="rad"))        ; t["phi"].axis_label = r"$\phi$"
    t.add_column(table.Column(data = as_astrometry[1], name="theta", unit="rad"))      ; t["theta"].axis_label = r"$\theta$"
    t.add_column(table.Column(data = as_astrometry[2], name="plx", unit="mas"))        ; t["plx"].axis_label = r"$\bar\omega$"
    t.add_column(table.Column(data = as_astrometry[3], name="muphi*", unit="mas/yr"))  ; t["muphi*"].axis_label = r"$\mu_\phi^*$"
    t.add_column(table.Column(data = as_astrometry[4], name="mutheta", unit="mas/yr")) ; t["mutheta"].axis_label = r"$\mu_\theta$"
    t.add_column(table.Column(data = as_astrometry[5], name="v_r", unit="km/s"))       ; t["v_r"].axis_label = "$v_r$"

    t["plx"] = t["plx"].to("arcsec")
    t["plx"].unit = gen.units.arcsec
    t["muphi*"] = t["muphi*"].to("arcsec/yr")
    t["muphi*"].unit = gen.units.arcsec / gen.units.year
    t["mutheta"] = t["mutheta"].to("arcsec/yr")
    t["mutheta"].unit = gen.units.arcsec / gen.units.year

def add_A(t):
    all_A = gen.map_np(gen.A, t["phi"], t["theta"])

    t.add_column(table.Column(data = all_A, name="A"))

def add_w(t):
    kplx = gen.k / t["plx"]
    w1 = np.zeros_like(t["v_r"])
    w2 = kplx * t["muphi*"] # muphi* = muphi * cos theta
    w3 = kplx * t["mutheta"]

    t.add_column(table.Column(data = w1, name = "w1", unit = "km / s")) ; t["w1"].axis_label = r"$w_1$"
    t.add_column(table.Column(data = w2, name = "w2", unit = "1 / yr")) ; t["w2"].axis_label = r"$w_2$"
    t.add_column(table.Column(data = w3, name = "w3", unit = "1 / yr")) ; t["w3"].axis_label = r"$w_3$"

def plot_distances(t, saveto):
    dist = gen.r(t["x"], t["y"], t["z"])
    plt.hist(dist, bins=25)
    plt.xlabel("Distance (pc)")
    plt.ylabel("Number of stars")
    plt.title("Distance distribution of simulated stars")
    plt.axis("tight")
    plt.tight_layout()
    plt.savefig(saveto)
    plt.close()

def plot_columns(t, saveto, *cols):
    cols_not_in_t = [col for col in cols if col not in t.keys()]
    if cols_not_in_t:
        raise gen.GFC_Exception("\n\ngaia_fc.simulated.plot_columns: The following columns are missing in the table: {0}".format(cols_not_in_t))
    for j, p1 in enumerate(cols):
        for p2 in cols[j+1:]:
            plt.scatter(t[p1], t[p2], c = t["group"], **p)
            plt.xlabel(t[p1].axis_label + " (" + t[p1].unit.to_string() + ")")
            plt.ylabel(t[p2].axis_label + " (" + t[p2].unit.to_string() + ")")
            plt.title("Simulated stars distribution: {0} vs {1}".format(p1, p2))
            plt.axis("tight")
            plt.tight_layout()
            plt.savefig(saveto.format(p1, p2))
            plt.close()

def plot_pdf_and_v_data(t, saveto, PDF, title="{0} - {1}"):
    "hard coded to work with v_x, v_y, v_z"

    cols_not_in_t = [col for col in ("v_x", "v_y", "v_z", "group") if col not in t.keys()]
    if cols_not_in_t:
        raise gen.GFC_Exception("\n\ngaia_fc.simulated.plot_pdf_and_v_data: The following columns are missing in the table: {0}".format(cols_not_in_t))
    if PDF.ndim != 3:
        raise gen.GFC_Exception("\n\ngaia_fc.simulated.plot_pdf_and_v_data: The provided PDF has dimensions {0} -- expected 3 dimensions".format(PDF.ndim))

    for i, pair in enumerate((("v_y", "v_z"), ("v_x", "v_z"), ("v_x", "v_y"))):
        p1, p2 = pair
        pair_int = [0, 1, 2]
        pair_int.remove(i)
        pdfsumi = PDF.sum(axis=i).T
        plt.scatter(t[p1], t[p2], c=t["group"], **p)
        plt.axis("tight")
        plt.contour(pdfsumi, extent=extent, linewidths=3, linestyles="dashed", levels=levels, colors='k')
        plt.xlabel(t[p1].axis_label + " (" + t[p1].unit.to_string() + ")")
        plt.ylabel(t[p2].axis_label + " (" + t[p2].unit.to_string() + ")")
        plt.title(title.format(p1, p2))
        plt.tight_layout()
        plt.savefig(saveto.format(p1, p2))
        plt.close()

def convert_back_from_astrometry_without_v_r(t):
    cols_already_in_t = [col for col in ("v_x_b", "v_y_b", "v_z_b") if col in t.keys()]
    if cols_already_in_t:
        raise gen.GFC_Exception("\n\ngaia_fc.simulated.convert_back_from_astrometry_without_v_r: The following columns are already in the table: {0}".format(cols_already_in_t))

    as_phase = np.array(gen.tophase(t["phi"], t["theta"], 1000.*t["plx"], 1000.*t["muphi*"], 1000.*t["mutheta"], np.zeros_like(t["v_r"])))
    t.add_column(table.Column(data = as_phase[3], name="v_x_B", unit="km/s")); t["v_x_B"].axis_label = r"$v_x (B)$"
    t.add_column(table.Column(data = as_phase[4], name="v_y_B", unit="km/s")); t["v_y_B"].axis_label = r"$v_y (B)$"
    t.add_column(table.Column(data = as_phase[5], name="v_z_B", unit="km/s")); t["v_z_B"].axis_label = r"$v_z (B)$"

def reconstruct_v_from_astrometry(t):
    cols_already_in_t = [col for col in ("v_x_R", "v_y_R", "v_z_R") if col in t.keys()]
    if cols_already_in_t:
        raise gen.GFC_Exception("\n\ngaia_fc.simulated.reconstruct_v_from_astrometry: The following columns are already in the table: {0}".format(cols_already_in_t))

    w = np.vstack((t["w1"], t["w2"], t["w3"])).T
    reconstr = gen.map_np(lambda A, w: A.dot(w), t["A"], w)

    t.add_column(table.Column(data = reconstr[:, 0], name = "v_x_R", unit = "km / s")); t["v_x_R"].axis_label = r"$v_x (R)$"
    t.add_column(table.Column(data = reconstr[:, 1], name = "v_y_R", unit = "km / s")); t["v_y_R"].axis_label = r"$v_y (R)$"
    t.add_column(table.Column(data = reconstr[:, 2], name = "v_z_R", unit = "km / s")); t["v_z_R"].axis_label = r"$v_z (R)$"
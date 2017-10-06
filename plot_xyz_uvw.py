"""
Plot all data.
Test - not used.
"""

import gfc
from gfc import ArgumentParser
from matplotlib import pyplot as plt
parser = ArgumentParser()
parser.add_argument("data_file", help = "File that contains the data")
parser.add_argument("save_folder", help = "Folder in which results will be saved")
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

t = gfc.io.read_csv(args.data_file)

gfc.add_rad(t)
x1,x2,x3,v1,v2,v3 = gfc.matrix.tophase(t["ra_rad"], t["dec_rad"], t["parallax"], t["pmra"], t["pmdec"], t["HRV"])
X, Y, Z = gfc.ICRS_to_galactic.transformCartesianCoordinates(x1, x2, x3)
U, V, W = gfc.ICRS_to_galactic.transformCartesianCoordinates(v1, v2, v3)
fig, axs = plt.subplots(1, 2)
im = axs[0].hexbin(-X, Z, W, gridsize = 40, vmin = -50, vmax = 50)
axs[1].hexbin(8000.-X, Z, W, gridsize = 40, vmin = -50, vmax = 50, extent = (5500, 10500, -2100, 2100))
plt.colorbar(im)
plt.show()

import gfc
import numpy as np
import copy as cp
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from astropy.table import Table

root = '/disks/strw9/vanweenen/mrp1/data/'

hipfile = 'Gaia_DR1_hipparcos_plx:eplx<0.1_nobinary'
tgasfile = 'Gaia_DR1_tgas_2mass_qualAAA_plx:eplx<0.1'

print "Reading files.."
print "Original.."
t_hip = gfc.io.read_csv(root + 'original/' + hipfile + '.csv')
t_tgas = gfc.io.read_csv(root + 'original/' + tgasfile + '.csv')

print "Main Sequence cut.."
t_tgas_mscut = gfc.io.read_csv(root + 'MS cut/' + tgasfile + '_MScut.csv')
t_tgas_mscut2 = gfc.io.read_csv(root + 'MS cut/' + tgasfile + '_MScut2.csv')
t_tgas_mscut3 = gfc.io.read_csv(root + 'MS cut/' + tgasfile + '_MScut3.csv' )
t_hip_copy = cp.copy(t_hip)
print "Done"

# Data selection TGAS
def density(x, y, mscut = False):
    saveto = '/disks/strw9/vanweenen/mrp1/results/data selection/'
    filename = '/CMD'
    title = 'TGAS sample'
    if not mscut:
        title += ' before main sequence cut'
        filename += '.png'
    else:
        title += ' final'
        filename += '_mscut.png'
    plt.figure(figsize=(5,5), tight_layout=True)
    plt.hist2d(x, y, bins=175, cmap=plt.cm.viridis, norm = LogNorm()) #histkwargs = {"vmin":1, "vmax":1200}
    plt.xlim((-0.5,4.0)) ; plt.ylim((11,-1.8))
    #plt.gca().invert_yaxis()
    plt.xlabel("$G - K$") ; plt.ylabel("$G$")
    plt.title(title)
    plt.savefig(saveto+filename)
    plt.show()

density(t_tgas['g_min_ks'], t_tgas['g_mag_abs'], mscut=False)
density(t_tgas_mscut['g_min_ks'], t_tgas_mscut['g_mag_abs'], mscut=True)
density(t_tgas_mscut2['g_min_ks'], t_tgas_mscut2['g_mag_abs'], mscut=True)
density(t_tgas_mscut3['g_min_ks'], t_tgas_mscut3['g_mag_abs'], mscut=True)


# Data selection for Hipparcos
print "Plot CMD.."
def AumerBinneyCut(xmin, xmax):
    xcut = np.linspace(xmin, xmax)
    yupper = np.empty(len(xcut))
    ylower = np.empty(len(xcut))
    for i, x in enumerate(xcut): 
        if x <= .5:
            yupper[i] = 7.5 * x - 3.75
        elif x >= .5 and x <= .8:
            yupper[i] = 15.33 * x - 7.665
        elif x >= .8:
            yupper[i] = 4.43 * x + 1.055
        if x <= .35:
            ylower[i] = 4.62 * x + 2.383
        elif x >= .35 and x <= .65:
            ylower[i] = 8.33 * x + 1.0845
        elif x >= .65 and x <= 1.25:
            ylower[i] = 3.33 * x + 4.3375
        elif x >= 1.25:
            ylower[i] = 6.50 * x + 0.375
    return xcut, yupper, ylower
    
def plot_Hipparcos(x, y, mscut = False):
    saveto = '/disks/strw9/vanweenen/mrp1/results/data selection/'
    filename = '/CMD_Hipparcos'
    title = 'Hipparcos sample'
    if not mscut:
        title += ' before main sequence cut'
        filename += '.png'
    plt.figure(figsize=(5,5), tight_layout=True)
    plt.scatter(x, y, marker='.', c='black', s=.05)
    if mscut:
        xcut, yupper, ylower = AumerBinneyCut(min(x), max(x))
        plt.plot(xcut, yupper, c='black')
        plt.plot(xcut, ylower, c='black')
        filename += '_mscut.png'
    plt.gca().invert_yaxis()
    plt.xlabel("$B - V$ (mag)") ; plt.ylabel("$H_p$ (mag)")
    plt.title(title)
    plt.savefig(saveto+filename)
    plt.show()

def AbsoluteMagnitude(m, plx):
    """
    Calculate absolute magnitude
    Input:
        m - apparent magnitude
        plx (arcsecond) - parallax
    """
    return m + 5*np.log10(plx) - 10

hp_mag_abs = AbsoluteMagnitude(t_hip['hp_mag'], t_hip['plx'])

plot_Hipparcos(t_hip['b_v'], hp_mag_abs, mscut=False)
plot_Hipparcos(t_hip['b_v'], hp_mag_abs, mscut=True)

print "Data selection.."

def Selection_Hipparcos_AumerBinney(t, hp_col='hp_mag', bv_col='b_v', plx_col='plx'):
    """
    Select main sequence stars of Hipparcos catalogue using the colour-magnitude diagram cuts from Aumer and Binney (2009)
    """
    hp = t[hp_col]
    bv = t[bv_col]
    plx = t[plx_col]
    hp_abs = AbsoluteMagnitude(hp, plx)
    print "1", len(t)
    t.remove_rows(np.where(np.logical_not((hp_abs < 7.50 * bv - 3.75) & (bv <= .5))))
    print "2", len(t)
    t.remove_rows(np.where(np.logical_not((hp_abs < 15.33 * bv - 7.665) & np.logical_and(bv >= .5, bv <= .8))))
    print "3", len(t)
    t.remove_rows(np.where(np.logical_not((hp_abs < 4.43 * bv + 1.055) & (bv >= .8))))
    print "4", len(t)
    t.remove_rows(np.where(np.logical_not((hp_abs > 4.62 * bv + 2.383) & (bv <= .35))))
    print "5", len(t)
    t.remove_rows(np.where(np.logical_not((hp_abs > 8.33 * bv + 1.0845) & np.logical_and(bv >= .35, bv <= .65))))
    print "6", len(t)
    t.remove_rows(np.where(np.logical_not((hp_abs > 3.33 * bv + 4.3375) & np.logical_and(bv >= .65, bv <= 1.25))))
    print "7", len(t)
    t.remove_rows(np.where(np.logical_not((hp_abs > 6.50 * bv + 0.375) & (bv >= 1.25))))
    return t

print len(t_hip_copy)
t_hip_copy = Selection_Hipparcos_AumerBinney(t_hip_copy)
print len(t_hip_copy)
t_hip_copy = cp.copy(t_hip)


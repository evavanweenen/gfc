"""
Functions to more easily make certain plots.
Some are very general, some are very hardcoded.
"""

from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec
viridis = plt.cm.viridis
from matplotlib import rc
rc('font',**{'family':'serif', 'serif':['Times']})
rc('text', usetex=True)

import numpy as np
import scipy as sp

from scipy.stats import chi2

from operator import add

from itertools import cycle
symbolwheel = cycle(['o', 'v', 's', '^'])

def show_or_save(saveto = None, fig = None, *args, **kwargs):
    if saveto is None: # show
        try:
            plt.show(fig, *args, **kwargs)
        except AttributeError:
            plt.show(*args, **kwargs)
    else:
        try:
            fig.savefig(saveto, *args, **kwargs)
        except AttributeError:
            plt.savefig(saveto, *args, **kwargs)

def title_filename(title, filename, *args):
    print("args:", args)
    if args:        
        if args[0] == 'Bovy':
            print("Bovy args detected")
            title += ' for Bovy (2009) input parameters'
            filename += '_Bovy'
        elif args[0] == 'Gaia':
            print("Gaia args detected")
            title += ' for Gaia data'
            filename += '_Gaia'
        else:
            param_title = ('K', 'N', '$\Delta\mu$', '$\sigma$', 'p', 'q')
            param_file = ('K', 'N', 'mean', 'sigma', 'p', 'q')
            title += ' with'
            for arg in range(len(args)):
                title += ' %s = %s' %(param_title[arg], args[arg])
                filename += '_%s%s' %(param_file[arg], args[arg])
        try:
            if args[1] == 'XD':
                title += ' fitted with XD'
                filename += '_logL'
            elif args[1] == 'AIC':
                title += ' fitted with AIC'
                filename += '_AIC'
            elif args[1] == 'MDL':
                title += ' fitted with MDL'
                filename += '_MDL'
        except IndexError:
            pass
    filename +='.png'
    return title, filename

def plot_normal_PDF(v, x, ax, amps, means, covs, c, l):
    pdf = np.zeros(len(x))
    for n in range(len(amps)):
        pdf += amps[n]*sp.stats.norm.pdf(x, loc=means[n,v], scale=np.sqrt(covs[n,v,v]))
    ax.plot(x, pdf, label=l, color=c, lw=.8)

def plot_hist_uvw(leny, initamps, initmeans, initcovs, a_test, m_test, c_test, uvw_data, a_test_vrad0=None, m_test_vrad0=None, c_test_vrad0=None, uvw_data_vrad0=None, *args):
    saveto = '/disks/strw9/vanweenen/mrp1/results/'
    velocities = ('U', 'V', 'W')
    unit = ' (km/s)'
    colors_vrad0 = ('red','dodgerblue','green')
    colors = ('lightcoral', 'skyblue', 'greenyellow')
    test = ('logL', 'AIC', 'MDL')
    test_vrad0 = ('logL $v_r = 0$', 'AIC $v_r = 0$', 'MDL $v_r = 0$')
    fig, ax = plt.subplots(int(leny), len(velocities), sharey=True, figsize=(9,3*int(leny)), tight_layout=True)
    def plots(ax, v, j):
        plot_normal_PDF(v, x, ax, initamps, initmeans, initcovs, 'black', 'initial')
        ax.hist(uvw_data[:,v], bins='auto', normed=True, facecolor='black', histtype='stepfilled', alpha=0.3, label='data')
        if uvw_data_vrad0 is not None:
            ax.hist(uvw_data_vrad0[:,v], bins='auto', normed=True, facecolor='grey', histtype='stepfilled', alpha=0.15, label='data $v_r = 0$')
        ax.set_xlabel(velocities[v] + unit)
        ax.legend(loc='upper right', prop={'size': 6})
    for v in range(len(velocities)):
        for j in range(len(test)):
            x = np.linspace(np.amin(uvw_data[:,v]), np.amax(uvw_data[:,v]), len(uvw_data[:,v]))
            if leny == 3:
                axx = ax[j,v]
                plot_normal_PDF(v, x, axx, a_test[j], m_test[j], c_test[j], colors[j], test[j])
                if a_test_vrad0 is not None:
                    plot_normal_PDF(v, x, axx, a_test_vrad0[j], m_test_vrad0[j], c_test_vrad0[j], colors_vrad0[j], test_vrad0[j])
                plots(axx, v, j)
            elif leny == 1:
                axx = ax[v]
                plot_normal_PDF(v, x, axx, a_test[j], m_test[j], c_test[j], colors[j], test[j])
                if a_test_vrad0 is not None:
                    plot_normal_PDF(v, x, axx, a_test_vrad0[j], m_test_vrad0[j], c_test_vrad0[j], colors_vrad0[j], test_vrad0[j])
        if leny == 1:
            plots(axx, v, j)
    suptitle = 'Histogram of velocity in Cartesian coordinates'
    filename = '/hist_velocity'
    suptitle, filename = title_filename(suptitle, filename, *args)
    plt.suptitle(suptitle, y=1., fontsize=12)
    plt.savefig(saveto + filename)
    plt.show()
         
def plot_XD_w_K(logL, AIC, MDL, bestK, bestw, vrad0, Kmin, Kmax, wmin, wmax, *args):
    saveto = '/disks/strw9/vanweenen/mrp1/results/'
    title = ("logL","AIC" ,"MDL")
    test = (logL, AIC, MDL)
    c = ("black", "white", "white")
    fig, ax = plt.subplots(1,len(test),figsize=(3*len(test)+1,3),tight_layout=True)
    for i in range(len(test)):
        ax[i].set_xlabel("$K$") ; ax[i].set_ylabel("$\sqrt{w}$ (km s$^{-1}$)")
        ax[i].set_xlim(Kmin-.5, Kmax+.5) ; ax[i].set_ylim(wmin-.125, wmax+.125)
        ax[i].set_title(title[i])
        ax[i].scatter(bestK[i], np.sqrt(bestw[i]), color=c[i], marker="o", s=150)
        im = ax[i].imshow(test[i].T, cmap=plt.cm.viridis, extent=(Kmin-.5, Kmax+.5, wmin-.125, wmax+.125), aspect='auto', interpolation='none', origin="lower")
        cbar = fig.colorbar(im, ax=ax[i], orientation="vertical")#, format='%.0e')
    suptitle = 'Extreme deconvolution'
    filename = '/logL-Kw'
    if vrad0:
        suptitle += ' for $v_r = 0$'
        filename += '_vrad0'
    suptitle, filename = title_filename(suptitle, filename, *args)
    plt.suptitle(suptitle, y=1., fontsize=11)
    plt.savefig(saveto + filename)
    plt.show()

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
    filename = '/CMD_Hipparcos.png'
    title = 'Hipparcos sample'
    if not mscut:
        title += ' before main sequence cut'
    plt.figure(figsize=(5,5), tight_layout=True)
    plt.scatter(x, y, c='black', s=.1)
    if mscut:
        xcut, yupper, ylower = AumerBinneyCut(min(x), max(x))
        plt.plot(xcut, yupper, c='black')
        plt.plot(xcut, ylower, c='black')
    plt.gca().invert_yaxis()
    plt.xlabel("$B - V$ (mag)") ; plt.ylabel("$H_p$ (mag)")
    plt.title(title)
    plt.savefig(saveto+filename)
    plt.show()

def density(x, y, mscut = False):
    saveto = '/disks/strw9/vanweenen/mrp1/results/data selection/'
    filename = '/CMD.png'
    title = 'TGAS sample'
    if not mscut:
        title += ' before main sequence cut'
    plt.figure(figsize=(5,5), tight_layout=True)
    plt.hist2d(x, y, bins=175, cmap=plt.cm.viridis, norm = LogNorm()) #histkwargs = {"vmin":1, "vmax":1200}
    plt.xlim((-0.5,4.0)) ; plt.ylim((11,-1.8))
    #plt.gca().invert_yaxis()
    plt.xlabel("$G - K$") ; plt.ylabel("$G$")
    plt.title(title)
    plt.savefig(saveto+filename)
    plt.show()

def eigsorted(cov):
    """
    http://www.nhsilbert.net/source/2014/06/bivariate-normal-ellipse-plotting-in-python/
    """
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]

def draw_PDF_ellipse(ax, amp, mean, cov, xyz, volume = 0.6827, **kwargs):
    choose_from = {"xy": [0, 1], "xz": [0, 2], "yz": [1, 2]}
    pair_int = choose_from[xyz]
    mean_ = mean[pair_int]
    cov_ = cov[pair_int][:, pair_int]

    vals, vecs = eigsorted(cov_)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * np.sqrt(chi2.ppf(volume,2)) * np.sqrt(vals)
    ell = Ellipse(xy = mean_, width = width, height = height, angle = theta, facecolor = "None", linewidth = 1 + 7.5*amp, **kwargs)
    ax.add_artist(ell)

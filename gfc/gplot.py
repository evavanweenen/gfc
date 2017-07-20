"""
Functions to more easily make certain plots.
Some are very general, some are very hardcoded.
"""

from matplotlib import pyplot as plt
from matplotlib.pyplot.cm import viridis
from matplotlib.patches import Ellipse
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec

import numpy as np

from scipy.stats import chi2

from .pdf import entropy

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

def PDFs(PDFs, axis, saveto = None, volume = 0.6827, ellipse_kwargs = {}, **kwargs):
    assert axis in ("xy", "xz", "yz"), "gaia_fc.gplot.PDFs: parameter `axis` must be one of ('xy', 'xz', 'yz'); received value {0}".format(axis)
    plt.figure()
    ax = plt.gca()
    for pdf in PDFs:
        draw_PDF_ellipse(pdf, ax, axis, volume = volume, **ellipse_kwargs)
    if kwargs:
        plt.setp(ax, **kwargs)
    plt.tight_layout()
    show_or_save(saveto)
    plt.close()

def PDFs_multi(PDFs, saveto = None, volume = 0.6827, title = "Gaussians", **kwargs):
    try:
        for i, pair, xl, yl in zip([0, 1, 2], [("v_y", "v_z"), ("v_x", "v_z"), ("v_x", "v_y")], [(-130, 120), (-130, 120), (-130, 120)], [(-70, 70), (-70, 70), (-120, 60)]):
            plt.figure()
            ax = plt.gca()
            for pdf in PDFs:
                p1, p2 = pair
                pair_int = [0, 1, 2]
                pair_int.remove(i)

                mean = pdf.mean[pair_int]
                cov = pdf.cov[pair_int][:, pair_int]

                vals, vecs = eigsorted(cov)
                theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
                width, height = 2 * np.sqrt(chi2.ppf(volume,2)) * np.sqrt(vals)
                ell = Ellipse(xy = mean, width = width, height = height, angle = theta, facecolor = "None", linewidth = 1 + 7.5*pdf.amp, **kwargs)
                ax.add_artist(ell)

            plt.xlabel("${0}$ (km/s)".format(pair[0]))
            plt.ylabel("${0}$ (km/s)".format(pair[1]))
            plt.title(title.format(p1, p2))
            plt.xlim(xl)
            plt.ylim(yl)
            if kwargs:
                plt.setp(**kwargs)
            plt.tight_layout()
            if saveto is None:
                plt.show()
            else:
                plt.savefig(saveto.format(p1, p2))
            plt.close()
    finally:
        plt.close("all")

def select_colourmap(name=None, default="jet"):
    try:
        cm = plt.cm.__dict__[name]
    except:
        cm = plt.cm.__dict__[default]
    return cm

def PDFs_gradients(PDFs_evaluated, saveto=None, log = True, cb = False, cmap="jet", extent = None, vmin=None, vmax=None, origin="lower", imshowkwargs={}, **kwargs):
    plt.figure()
    cm = select_colourmap(cmap)
    norm = LogNorm() if log else None
    plt.imshow(PDFs_evaluated.T, cmap=cm, extent = extent, vmin=vmin, vmax=vmax, norm=norm, origin=origin, **imshowkwargs)
    plt.grid(False)
    plt.axis("tight")
    if kwargs:
        plt.setp(plt.gca(), **kwargs)
    plt.tight_layout()
    show_or_save(saveto)
    plt.close()

def density(x, y, saveto=None, bins=500, cb = False, flip="", log = True, cmap="viridis", r=None, grid=False, histkwargs={}, **kwargs):
    norm = LogNorm() if log else None
    cm = select_colourmap(cmap)
    plt.figure()
    try:
        _ = plt.hist2d(x, y, bins=bins, cmap=cm, range = r, norm = norm, **histkwargs) ; del _
        if cb:
            plt.colorbar(pad = 0.01)
        plt.grid(grid)
        plt.axis("tight")
        if "x" in flip:
            plt.gca().invert_xaxis()
        if "y" in flip:
            plt.gca().invert_yaxis()
        try:
            plt.xlabel(x.name)
            plt.ylabel(y.name)
        except:
            pass
        if kwargs:
            plt.setp(plt.gca(), **kwargs)
        plt.tight_layout()
        show_or_save(saveto)
    finally:
        plt.close()

def density_ax(ax, x, y, bins=500, cb = False, flip="", log = True, cmap="viridis", r=None, grid=False, histkwargs={}, **kwargs):
    norm = LogNorm() if log else None
    cm = select_colourmap(cmap)
    _ = ax.hist2d(x, y, bins=bins, cmap=cm, range = r, norm = norm, **histkwargs) ; del _
    #if cb:
    #    plt.colorbar(pad = 0.01)
    ax.grid(grid)
    ax.axis("tight")
    if "x" in flip:
        ax.invert_xaxis()
    if "y" in flip:
        ax.invert_yaxis()
    try:
        ax.set_xlabel(x.name)
        ax.set_ylabel(y.name)
    except:
        pass
    if kwargs:
        plt.setp(ax, **kwargs)

def RV_distribution(v_r, p_v_r, saveto=None, name=None, real_vr=None, entropy_text=True, labels=None, **kwargs):
    try:
        plt.figure()
        if len(p_v_r.shape) == 1: # only one distribution
            plt.plot(v_r, p_v_r, c='k', lw=2, label="$P(v_r)$")
        else:
            if labels is None:
                labels = ["${0}$".format(k) for k, p in enumerate(p_v_r)]
            for p, l in zip(p_v_r, labels):
                plt.plot(v_r, p, lw=2, label=l)
        if real_vr is not None:
            plt.axvline(real_vr, c='k', lw=2, ls="--", label="Observed")
        plt.xlabel(r"$v_r$")
        plt.ylabel(r"$P (v_r)$")
        plt.axis("tight")
        if name is not None:
            plt.text(0.05, 0.95, str(name), transform=plt.gca().transAxes, ha="left")
        if entropy_text:
            if len(p_v_r.shape) == 1:
                entropies = ["$H_1$: {0:.3f}".format(entropy(v_r, p_v_r))]
            else:
                entropies = ["$H_{0}$: {1:.3f}\n".format(k+1, entropy(v_r, p)) for k, p in enumerate(p_v_r)]
            entropiestr = reduce(add, entropies)
            plt.text(0.7, 0.95, entropiestr, transform=plt.gca().transAxes, va="top")
        plt.legend(loc="center right")
        if kwargs:
            plt.setp(plt.gca(), **kwargs)
        plt.tight_layout()
        show_or_save(saveto)
    finally:
        plt.close()
        
def plotvr(vr, pvr, saveto=None, *args, **kwargs):
    try:
        plt.figure(figsize=(10,10), tight_layout=True)

        plotvrax(plt.gca(), vr, pvr, *args, **kwargs)
                
        show_or_save(saveto)
    finally:
        plt.close()
    
def plotvrax(ax, vr, pvr, name=None, HRV=None, eHRV=None, lower_c=None, upper_c=None, H_c=None, H_m=None, **kwargs):
    ax.plot(vr, pvr[0], lw=2, c='r', label = '$p(v_r|\mathbf{v}_t)$')
    ax.plot(vr, pvr[1], lw=2, c='b', label = '$p(v_r)$')

    if HRV is not None:
        ax.axvline(HRV, label = '$v_r$ (RAVE)', c='k', lw=2, ls='--')
        if eHRV is not None:
            ax.axvspan(HRV-eHRV, HRV+eHRV, color='0.5', alpha=0.5)
    if lower_c is not None and upper_c is not None:
        ax.axvline(lower_c, c='r', ls='--', lw=2)
        ax.axvline(upper_c, c='r', ls='--', lw=2)

    ax.axis("tight")
    ax.set_ylim(0, pvr.max()*1.05)
    ax.set_xlim(vr[0]-1, vr[-1]+1)

    if H_c is not None and H_m is not None:
        entropiestr = "$H_C = {hc:.3f}$\n$H_M = {hm:.3f}$".format(hc = H_c, hm = H_m)
        ax.text(0.9, 0.9, entropiestr, transform=ax.transAxes, va="top", ha="right")
    if name is not None:
        ax.text(0.05, 0.9, name, transform=ax.transAxes, ha="left")
        
    if kwargs:
        ax.set(**kwargs)
    
def binned_1d_x_xerr(bin_left_edges, highest):
    left = np.concatenate((np.array(bin_left_edges), np.array([highest])))
    center = np.array([(this+previous)/2. for this, previous in zip(left[1:], left[:-1])])
    xerr = center - left[:-1]
    return center, xerr

def binned_1d(bin_left_edges, values, highest, sigmas = None, saveto=None, x_sig=None, dotkwargs={}, **kwargs):
    center, xerr = binned_1d_x_xerr(bin_left_edges, highest)
    if x_sig is not None:
        xerr = x_sig
    plt.figure()
    plt.errorbar(center, values, xerr=xerr, yerr=sigmas, fmt="o", c='k', **dotkwargs)
    plt.axis("tight")
    if kwargs:
        plt.setp(plt.gca(), **kwargs)
    show_or_save(saveto)
    plt.close()

def binned_1d_multi(bin_left_edges_list, values_list, highest_list, sigmas_list=None, saveto=None, kwargs_list=None, dotkwargs={}, fmt="o", x_sig=None, **figure_kwargs):
    """
    this function is a MESS
    should be redone completely
    if x_sig then sigma_list is assumed to be zip(xerr, yerr)
    """
    fig, axs = plt.subplots(ncols=1, nrows=len(bin_left_edges_list), **figure_kwargs)
    try:
        if sigmas_list is None:
            sigmas_list = [None]*len(bin_left_edges_list)
        if kwargs_list is None:
            kwargs_list = [None]*len(bin_left_edges_list)
        if x_sig is None:
            x_sig = [None]*len(bin_left_edges_list)
        assert len(bin_left_edges_list) == len(values_list) == len(highest_list) == len(sigmas_list) == len(kwargs_list), "gaia_fc.gplot.binned_1d_multi: Not all inputs are of the same length: {0}-{1}-{2}-{3}-{4}".format([len(x) for x in (bin_left_edges_list, values_list, highest_list, sigmas_list, kwargs_list)])
        for a,b,v,h,s,kw,x_ in zip(axs, bin_left_edges_list, values_list, highest_list, sigmas_list, kwargs_list, x_sig):
            if any(isinstance(Z, tuple) for Z in (b, v, h, s, x_)): # multiple scatters in one axes subplot
                assert all(isinstance(Z, tuple) for Z in (b, v, h, s)), "gaia_fc.gplot.binned_1d_multi: Not all objects to be plotted are of the same type: {0} {1} {2} {3}".format(*[type(Z) for Z in (b, v, h, s)])
                if not isinstance(x_, tuple):
                    x_ = [x_ for i in b]
                for b_, v_, h_, s_, x__, f in zip(b,v,h,s,x_,symbolwheel):
                    center, xerr = binned_1d_x_xerr(b_, h_)
                    if x__ is not None:
                        xerr = x__
                    a.errorbar(center, v_, xerr=xerr, yerr=s_, fmt=f, c='k', **dotkwargs)
            else:
                center, xerr = binned_1d_x_xerr(b, h)
                if x_ is not None:
                    xerr = x_
                a.errorbar(center, v, xerr=xerr, yerr=s, fmt=fmt, c='k', **dotkwargs)
            a.axis("tight")
            if kw:
                plt.setp(a, **kw)
        fig.tight_layout()
        show_or_save(saveto, fig)
    finally:
        plt.close(fig)
        
def BV_SUVW(UR, VR, WR, S2R, UM, VM, WM, S2M, saveto=None):
    # hard coded version of binned_1d_multi for ease of use
    f, axs = plt.subplots(3, 2, sharex=True, sharey="row", gridspec_kw={"height_ratios": (2,1,1)}, tight_layout = True)
    axsR = axs[:,0] ; axsM = axs[:,1]
    for axs, U, V, W, S2 in zip((axsR, axsM), (UR, UM), (VR, VM), (WR, WM), (S2R, S2M)):
        c, xerr = binned_1d_x_xerr(U[0], U[3])
        axs[1].errorbar(c, U[1], xerr=xerr, yerr=U[2], fmt="o",c="k")
        axs[0].errorbar(c, S2[1], xerr=xerr, yerr=S2[2], fmt="s", c="k", label = "$S$")
        axs[0].errorbar(c, V[1], xerr=xerr, yerr=V[2], fmt="^",c="k", label = "$V$")
        axs[0].legend(loc = "lower right", ncol=2, numpoints=1)
        axs[2].errorbar(c, W[1], xerr=xerr, yerr=W[2], fmt="v",c="k")
        axs[2].set_xlabel("$B - V$")
        axs[2].set_xlim(-0.3, 1.7)
        axs[2].set_xticks((0, 0.5, 1, 1.5))
    axsR[0].set_ylabel("$V, S$ (km s$^{-1})$")
    axsR[1].set_ylabel("$U$ (km s$^{-1})$")
    axsR[2].set_ylabel("$W$ (km s$^{-1})$")
    axsR[0].set_title("$v_r$ from RAVE")
    axsM[0].set_title("$v_r$ from XD")
    axsR[0].set_ylim(0, 45)
    axsR[1].set_ylim(0, 20)
    axsR[2].set_ylim(0, 20)
    show_or_save(saveto, f)
    plt.close(f)
    
def BV_SV(VR, S2R, VM, S2M, saveto=None):
    # hard coded version of binned_1d_multi for ease of use
    f, axs = plt.subplots(1, 2, sharex=True, sharey="row", tight_layout = True, figsize=(12, 7))
    axR = axs[0] ; axM = axs[1]
    for ax, V, S2 in zip(axs, (VR, VM), (S2R, S2M)):
        c, xerr = binned_1d_x_xerr(V[0], V[3])
        ax.errorbar(c, S2[1], xerr=xerr, yerr=S2[2], fmt="s", c="k", label = "$S$")
        ax.errorbar(c, V[1], xerr=xerr, yerr=V[2], fmt="^",c="k", label = "$V$")
        ax.legend(loc = "lower right", ncol=2, numpoints=1)
        ax.set_xlabel("$B - V$")
    axR.set_ylabel("$V, S$ (km s$^{-1})$")
    axR.set_title("$v_r$ from RAVE")
    axM.set_title("$v_r$ from XD")
    axR.set_ylim(0, 45)
    axR.set_xlim(-0.5, 2.3)
    axR.set_xticks((0, 0.5, 1, 1.5, 2))
    show_or_save(saveto, f)
    plt.close(f)

def S2UVW(UR, VR, WR, S2R, UM, VM, WM, S2M, UmeanR, VcoevR, WmeanR, UmeanM, VcoevM, WmeanM, saveto=None):
    f, axs = plt.subplots(3, 2, sharex=True, sharey="row", gridspec_kw={"height_ratios": (2,1,1)}, tight_layout = True)
    axsR = axs[:,0] ; axsM = axs[:,1]
    for axs, U, V, W, S2, Um, Vl, Wm in zip((axsR, axsM), (UR, UM), (VR, VM), (WR, WM), (S2R, S2M), (UmeanR, UmeanM), (VcoevR, VcoevM), (WmeanR, WmeanM)):
        axs[1].errorbar(S2[1]**2., U[1], yerr=U[2], fmt="o", c="k")
        axs[0].errorbar(S2[1]**2., V[1], yerr=V[2], fmt="^", c="k")
        axs[2].errorbar(S2[1]**2., W[1], yerr=W[2], fmt="v", c="k")
        axs[1].axhline(Um, ls="--", lw=2, c="k")
        axs[2].axhline(Wm, ls="--", lw=2, c="k")
        lx = np.arange(0, 2500, 500)
        ly = Vl[1] + lx * Vl[0]
        axs[0].plot(lx, ly, ls="--", lw=2, c="k")
        axs[2].set_xlabel("$S^2$ (km$^2$ s$^{-2}$)")
    axsR[0].set_ylabel("$V$ (km s$^{-1})$")
    axsR[1].set_ylabel("$U$ (km s$^{-1})$")
    axsR[2].set_ylabel("$W$ (km s$^{-1})$")
    axsR[0].set_title("$v_r$ from RAVE")
    axsM[0].set_title("$v_r$ from XD")
    axsR[0].set_ylim(0, 35)
    axsR[1].set_ylim(0, 20)
    axsR[2].set_ylim(0, 20)
    axsR[2].set_xlim(0, 2000)
    axsR[2].set_xticks((250, 500, 750, 1000, 1250, 1500, 1750))
    axsR[0].set_yticks((5,10,15,20,25,30))
    axsR[1].set_yticks((5,10,15))
    axsR[2].set_yticks((5,10,15))
    show_or_save(saveto, f)
    plt.close(f)

def S2V(VR, S2R, VM, S2M, VcoevR, VcoevM, saveto=None):
    # hard coded version of binned_1d_multi for ease of use
    f, axs = plt.subplots(1, 2, sharex=True, sharey="row", tight_layout = True, figsize=(12, 7))
    axR = axs[0] ; axM = axs[1]
    for ax, V, S2, Vl in zip(axs, (VR, VM), (S2R, S2M), (VcoevR, VcoevM)):
        ax.errorbar(S2[1]**2., V[1], yerr=V[2], fmt="^", c="k")
        lx = np.arange(0, 2500, 500)
        ly = Vl[1] + lx * Vl[0]
        ax.plot(lx, ly, ls="--", lw=2, c="k")
        ax.set_xlabel("$S^2$ (km$^2$ s$^{-2}$)")
    axR.set_ylabel("$V$ (km s$^{-1})$")
    axR.set_title("$v_r$ from RAVE")
    axM.set_title("$v_r$ from XD")
    axR.set_ylim(0, 35)
    axR.set_xlim(0, 2000)
    axR.set_xticks((250, 750, 1250, 1750))
    axR.set_yticks((5,10,15,20,25,30))
    show_or_save(saveto, f)
    plt.close(f)

def moving_group(ra, dec, l, b, U, V, W, ind_in, ind_ex = None, met = None, alpha = None, a = None, m = None, c = None, saveto = None, **kwargs):
    fig, axs = plt.subplots(2, 3, figsize = (15, 10), tight_layout = True)
    density_ax(axs[0,0], U[ind_in], V[ind_in], r = ((-130, 130), (-130, 130)), **kwargs)
    density_ax(axs[0,1], U[ind_in], W[ind_in], r = ((-130, 130), (-130, 130)), **kwargs)
    density_ax(axs[0,2], V[ind_in], W[ind_in], r = ((-130, 130), (-130, 130)), **kwargs)
    if a is not None and m is not None and c is not None:
        draw_PDF_ellipse(axs[0,0], a, m, c, "xy", zorder = 1000)
        draw_PDF_ellipse(axs[0,1], a, m, c, "xz", zorder = 1000)
        draw_PDF_ellipse(axs[0,2], a, m, c, "yz", zorder = 1000)
    density_ax(axs[1,0], ra[ind_in], dec[ind_in], r = ((0, 360), (-90, 90)),**kwargs)
    density_ax(axs[1,1], l[ind_in], b[ind_in], r = ((0, 360), (-90, 90)), **kwargs)
    if met is not None and alpha is not None:
        density_ax(axs[1,2], met[ind_in], alpha[ind_in], r = ((-2, 1), (-1, 1)), **kwargs)
    elif ind_ex is not None:
        density_ax(axs[1,2], l[ind_ex], b[ind_ex], r = ((0, 360), (-90, 90)), **kwargs)
    show_or_save(saveto, fig)

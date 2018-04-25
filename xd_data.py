import gfc
from gfc import *
import numpy as np
import scipy as sp
from scipy.stats import chi2
import copy as cp
from extreme_deconvolution import extreme_deconvolution as xd
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import gridspec as gs
from matplotlib.patches import Ellipse
import matplotlib as mpl
from mcmcplotting import convert_to_stdev_nan
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import Angle
from pygaia.astrometry.vectorastrometry import phaseSpaceToAstrometry, astrometryToPhaseSpace, sphericalToCartesian, normalTriad
from pygaia.astrometry.coordinates import CoordinateTransformation, Transformations
from pygaia.astrometry import constants

galtoicrs = CoordinateTransformation(Transformations.GAL2ICRS)
icrstogal = CoordinateTransformation(Transformations.ICRS2GAL)

rc('font',**{'family':'serif', 'serif':['Times']})
rc('text', usetex=True)

root = '/disks/strw9/vanweenen/mrp1/data/'

#hipfile = 'Gaia_DR1_hipparcos_plx:eplx<0.1_nobinary'
tgasfile = 'Gaia_DR1_tgas_2mass_qualAAA_plx:eplx<0.1_MScut'

print("Reading files..")
#t_hip = gfc.io.read_csv(root + 'original/' + hipfile + '.csv')
t_tgas = gfc.io.read_csv(root + 'MS cut/' + tgasfile + '.csv')

N = len(t_tgas)
print(N)

######################################################################################
# Use parameters from Bovy 2009 as initial velocity distribution.
amps_Bovy = gfc.io.load('/disks/strw9/vanweenen/mrp1/gfc/Bovy_parameters/' + 'Bovy_amps.npy')
means_Bovy = gfc.io.load('/disks/strw9/vanweenen/mrp1/gfc/Bovy_parameters/' + 'Bovy_means.npy')
covs_Bovy = gfc.io.load('/disks/strw9/vanweenen/mrp1/gfc/Bovy_parameters/' + 'Bovy_covs.npy')

amps_Bovy /= np.sum(amps_Bovy)

initamps = []
initmeans = []
initcovs = []

for k in range(len(amps_Bovy)):
    initamps.append(amps_Bovy[:k+1])
    initmeans.append(means_Bovy[:k+1])
    initcovs.append(covs_Bovy[:k+1])

Kmin=1 ; Kmax=20 ; Kstep=1
wmin = 0.5 ; wmax=15 ; wstep=1.
    
######################################################################################
print("Projection..")
# Calculate the projection matrix analogous to equation (1) in Bovy et al 2009 
# (https://ui.adsabs.harvard.edu/#abs/2009ApJ...700.1794B/abstract). Note the different ordering of the
# velocity components.
gfc.add_rad(t_tgas, 'ra', u.degree, errunit=u.mas)
gfc.add_rad(t_tgas, 'dec', u.degree, errunit=u.mas)

matrix.transformation(t_tgas)
warr = gfc.XD_arr(t_tgas, "w")
wcovar = gfc.XD_arr(t_tgas, "S") ; wcovar[:,2,2] = 1e15
proj = gfc.XD_arr(t_tgas, "R")
varr = gfc.XD_arr(t_tgas, "UVW")

print("Done")

# Perform XD
print("Perform XD..")
# The input to XD are the values of v_alpha*, v_delta, vrad. The other input required is the projection matrix.
# The values of v_alpha*, v_delta, vrad are obtained from (alpha, delta, parallax, mura,  mudec, vrad).
wrange = np.arange(wmin, wmax + wstep, wstep)**2.
Krange = range(Kmin, Kmax + Kstep, Kstep) #Krange from Kmin to Kmax

"""
logL, AIC, MDL, amps_test, means_test, covs_test, bestK, bestw = gfc.perform_XD(warr, wcovar, proj, initamps, initmeans, initcovs, wrange, Krange, N)

gfc.io.save_PDFs(amps_test, means_test, covs_test, root + '/Gaia')
"""
amps_test, means_test, covs_test = gfc.io.load_PDFs(root + '/Gaia/')

def title_filename(title, filename, *args, **kwargs):
    print("args:", args, "kwargs:", kwargs)
    if args:        
        for i in args:
            filename += '_'+i
        if args[0] == 'Bovy':
            print("Bovy args detected")
            title += ' for Bovy (2009) input parameters'
        elif args[0] == 'Gaia':
            print("Gaia args detected")
            title += ' for Gaia data'
        try:
            if args[1] == 'XD':
                title += ' fitted with XD'
            elif args[1] == 'AIC':
                title += ' fitted with AIC'
            elif args[1] == 'MDL':
                title += ' fitted with MDL'
        except IndexError:
            pass
    if kwargs:
        title += ' with'
        for i in kwargs:
            if i == 'meanstep':
                title += ' $\Delta\mu$ = %s km/s' %(kwargs[i])
            elif i == 'sigma':
                title += ' $\sigma$ = %s km/s' %(kwargs[i])
            else:
                title += ' %s = %s' %(i, kwargs[i])
            filename += '_%s%s' %(i, kwargs[i])
    filename +='.pdf'
    return title, filename

def plot_normal_PDF(ax, v, x, amps, means, covs, c, l):
    pdf = np.zeros(len(x))
    for n in range(len(amps)):
        pdf += amps[n]*sp.stats.norm.pdf(x, loc=means[n,v], scale=np.sqrt(covs[n,v,v]))
    ax.plot(x, pdf, label=l, color=c, lw=.8)

def plot_hist_uvw(inita, initm, initc, a_test, m_test, c_test, uvw_data, *args, **kwargs):
    saveto = '/disks/strw9/vanweenen/mrp1/results/'
    velocities = ('U', 'V', 'W') ; unit = ' (km/s)'
    colors = ('red','dodgerblue','green')
    test = ('logL', 'AIC', 'MDL')
    limits = ((-130, 120), (-120,60), (-70,70))
    
    fig, ax = plt.subplots(1, len(velocities), sharey=True, figsize=(10,4))
    fig.subplots_adjust(wspace=0.)
    for v in range(len(velocities)):
        vrange = np.linspace(limits[v][0], limits[v][1], len(uvw_data[:,v]))
        for t in range(len(test)):
            plot_normal_PDF(ax[v], v, vrange, a_test[t], m_test[t], c_test[t], colors[t], test[t]) #xd fit
        plot_normal_PDF(ax[v], v, vrange, inita, initm, initc, 'black', 'Bovy') #curve of Bovy initial values
        #ax[v].hist(uvw_data[:,v], bins='auto', normed=True, linewidth=2., edgecolor='black', histtype='step', alpha=0.3, label='data')
        ax[v].set_xlabel(velocities[v] + unit)
        ax[v].set_xlim(limits[v])
    ax[0].set_ylabel('Probability')
    ax[1].legend(loc='upper center', bbox_to_anchor=(0.5,1.10), ncol=5, fancybox=True, shadow=True, prop={'size': 8})
    suptitle = 'Histogram of velocity in Cartesian coordinates'
    filename = 'hist_velocity'
    suptitle, filename = title_filename(suptitle, filename, *args, **kwargs)
    plt.suptitle(suptitle, y=1.2, fontsize=12)
    plt.savefig(saveto + filename)
    plt.show()

plot_hist_uvw(amps_Bovy[:10], means_Bovy[:10], covs_Bovy[:10], amps_test, means_test, covs_test, varr, 'Gaia', N=N)

"""
def plot_XD_w_K(logL, AIC, MDL, bestK, bestw, vrad0, Kmin, Kmax, wmin, wmax, *args, **kwargs):
    saveto = '/disks/strw9/vanweenen/mrp1/results/'
    title = ("logL","AIC" ,"MDL")
    test = (logL, AIC, MDL)
    c = ("black", "white", "white")
    fig, ax = plt.subplots(1,len(test),figsize=(3*len(test)+1,3), sharey=True, tight_layout=True)
    for i in range(len(test)):
        ax[i].set_xlabel("$K$") 
        ax[i].set_xlim(Kmin-.5, Kmax+.5) ; ax[i].set_ylim(wmin-wstep/2., wmax+wstep/2.)
        ax[i].set_title(title[i])
        ax[i].scatter(bestK[i], np.sqrt(bestw[i]), color=c[i], marker="x", s=150)
        im = ax[i].imshow(test[i].T, cmap=plt.cm.viridis, extent=(Kmin-Kstep/2., Kmax+Kstep/2., wmin-wstep/2., wmax+wstep/2.), aspect='auto', interpolation='none', origin="lower")
        cbar = fig.colorbar(im, ax=ax[i], orientation="vertical")#, format='%.0e')
    ax[0].set_ylabel("$\sqrt{w}$ (km s$^{-1}$)")
    suptitle = 'Extreme deconvolution'
    filename = '/logL-Kw'
    if vrad0:
        suptitle += ' for $v_r = 0$'
        filename += '_vrad0'
    suptitle, filename = title_filename(suptitle, filename, *args, **kwargs)
    plt.suptitle(suptitle, y=0.999, fontsize=11)
    plt.savefig(saveto + filename)
    plt.show()

plot_XD_w_K(logL, AIC, MDL, bestK, bestw, False, Kmin, Kmax, wmin, wmax, 'Gaia')
"""

def KnownGroups(ax): #before: text
    ax.text(60, -90, "Arcturus")
    ax.plot([8, 57.5], [-105, -90], c='k', lw=2)
    ax.text(-80, -110, "Halo")
    ax.text(50, 40, "Sirius/UMa")
    ax.plot([49, 5], [40, 5], c='k', lw=2)
    ax.text(-100, 45, "Coma Berenices")
    ax.plot([-70, -10], [42, -5], c='k', lw=2)
    ax.text(-120, 34, "NGC 1901")
    ax.plot([-100, -25], [31, -12], c='k', lw=2)
    ax.text(-120, 0, "Hyades")
    ax.plot([-110, -45], [-3, -17], c='k', lw=2)
    ax.text(90, -50, "Pleiades")
    ax.plot([87, -15], [-45, -20], c='k', lw=2)
    ax.text(-125, -42, "Hercules")
    ax.plot([-93.5, -28], [-40, -42], c='k', lw=2)

def totalGaussians(axs, i, j, lsp, amps, means, covs):
    levels=np.array([2,6,12,21,33,50,68,80,90,95,99])/100.0
    xx, yy = np.meshgrid(lsp[j],lsp[i+1])
    pdfxy = np.zeros_like(xx)
    m = np.array([means[:,j], means[:,i+1]]).T
    c = np.array([[covs[:,j,j],covs[:,j,i+1]], [covs[:,i+1,j],covs[:,i+1,i+1]]]).T
    for k in range(len(amps)):
        pdfxy = pdfxy + amps[k] * \
        sp.stats.multivariate_normal.pdf(np.dstack((xx,yy)), mean=m[k], cov=c[k])
    axs[i+j].contour(lsp[j], lsp[i+1], convert_to_stdev_nan(np.log(pdfxy)), levels=levels, colors='k', linewidths=1)
                
def set_axes_3velocities(ax, limits):  
    velocities = ('U', 'V', 'W') ; unit = ' (km/s)'
    xpos = ('top', 'bottom') ; ypos = ('left', 'right')
    for i in range(2):
        for j in range(2):
            if j != i + 1:
                ax[i+j].set_xlim(limits[j]) ; ax[i+j].set_ylim(limits[i+1])
                ax[i+j].set_xlabel(velocities[j] + unit) ; ax[i+j].set_ylabel(velocities[i+1] + unit)
                ax[i+j].xaxis.set_label_position(xpos[i]) ; ax[i+j].yaxis.set_label_position(ypos[j])
    ax[0].xaxis.tick_top()
    ax[0].set_xticks((-100, -50, 0, 50, 100))
    ax[0].set_yticks((-100, -50, 0, 50))
    ax[2].yaxis.tick_right()

def plot_Gaussian_comps(amps, means, covs, t_uvw, edgecolor, total=False, showdata=False, showbovy=True, *args, **kwargs):
    saveto = '/disks/strw9/vanweenen/mrp1/results/'
    filename = '/PDF'
    limits = ((-130, 120), (-120,60), (-70,70)); gs = (251, 181, 161)
    if total:
        filename += '_total'
        title = "Velocity distribution"
        #levels = np.logspace(-6.,-2.7,10)
        lsp = [np.linspace(limits[i][0], limits[i][1], gs[i]) for i in range(3)]
    else:
        line = ("xy", "xz", "yz") 
        title = "Gaussian components of velocity distribution"
    title, filename = title_filename(title, filename, *args, **kwargs)
    fig = plt.figure(figsize=(8,8))
    fig.subplots_adjust(hspace=0.05) ; fig.subplots_adjust(wspace=0.05)
    plt.suptitle(title)
    ax = [plt.subplot(2,1,1), plt.subplot(2,2,3), plt.subplot(2,2,4)]
    for i in range(2):
        for j in range(2):
            if j != i +1:
                if total:
                    totalGaussians(ax, i, j, lsp, amps, means, covs)
                else:
                    if showbovy:
                        for a, m, c in zip(amps_Bovy[:10], means_Bovy[:10], covs_Bovy[:10]):
                            gfc.gplot.draw_PDF_ellipse(ax[i+j], a, m, c, line[i+j], edgecolor="0.4")
                    for a, m, c in zip(amps, means, covs):
                        gfc.gplot.draw_PDF_ellipse(ax[i+j], a, m, c, line[i+j], edgecolor=edgecolor)
                if showdata:
                    ax[i+j].hexbin(t_uvw[:,j], t_uvw[:,i+1], gridsize=gs[i+j], bins='log', mincnt=1, alpha=.4, lw=0.)
    set_axes_3velocities(ax, limits)    
    KnownGroups(ax[0])
    fig.savefig(saveto + filename)
    plt.show()
    
#plot_Gaussian_comps(amps_Bovy[:10], means_Bovy[:10], covs_Bovy[:10], varr, "0.4", False, True, False, 'Gaia', N=N)
#plot_Gaussian_comps(amps_Bovy[:10], means_Bovy[:10], covs_Bovy[:10], varr, "0.4", True, True, False, 'Gaia', N=N)
tests = ('XD', 'AIC', 'MDL')
edgecolors = ('red', 'blue', 'green')
for i in range(3):
    plot_Gaussian_comps(amps_test[i], means_test[i], covs_test[i], varr, edgecolors[i], False, False, True, 'Gaia', tests[i], N=N)
    plot_Gaussian_comps(amps_test[i], means_test[i], covs_test[i], varr, "0.4", True, False, False, 'Gaia', tests[i], N=N)

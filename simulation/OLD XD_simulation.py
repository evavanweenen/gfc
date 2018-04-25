
# coding: utf-8

# In[6]:

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


# In[7]:

root_folder = "/disks/strw9/vanweenen/mrp1/results/simulation"

# Simulate a collection of stars distributed at constant space density around the sun, 
# with distances between 1 and 500 pc. Use Galactic coordinates. Each star is assigned
# a 3D space velocity drawn from a superposition of several 3D normal distributions in 
# velocity space.
rmin = 1. #pc
rmax = 500. #pc
N = 1000 #nr of stars

r = np.random.uniform(rmin**3, rmax**3, N)**(1./3.) #pc
theta = np.arccos(np.random.uniform(-1, 1, N))-np.pi/2
phi = np.random.uniform(0, 2*np.pi, N)

x_gal,y_gal,z_gal = sphericalToCartesian(r, phi, theta)


# In[8]:

# Define the space velocities (u,v,w), again in galactic coordinates
K = 2 #nr of components
meanstep = 50. #km/s
sigma = 10. #km/s
p = 0. #relation between amplitudes of components
q = 0. #relation between sigmas of components

wmin = 0.5 ; wmax=30 ; wstep=1.5
#wmin = 0.5 ; wmax=5 ; wstep=.25
Kmin=1 ; Kmax=10 ; Kstep=1

def simulate_amps(k):
    a = np.empty([k])
    rest = 1
    for i in range(k-1):
        a[i] = (1./k)/(i+2)**p
        rest -= a[i]
    a[k-1] = rest
    return a

def simulate_means(k):
    m = np.empty([k, 3])
    for i in range(k):
        m[i,:] = [meanstep*i, meanstep*i, meanstep*i]
    return m

def simulate_covs(k):
    c = np.zeros((k,3,3))
    for i in range(k):
        c[i,0,0] = (sigma/((i+1)**q))**2
        c[i,1,1] = (sigma/((i+1)**q))**2
        c[i,2,2] = (sigma/((i+1)**q))**2
    return c

initamps = []
initmeans = []
initcovs = []

for k in range(Kmin, Kmax+1):
    initamps.append(simulate_amps(k))
    initmeans.append(simulate_means(k))
    initcovs.append(simulate_covs(k))

print initamps[K-1], initmeans[K-1], initcovs[K-1]

# Simulate velocity for every star
component = np.random.choice(K, N, p=initamps[K-1]) #choose component
uvw_gal = np.empty([N,3])
for i in range(N):
    uvw_gal[i] = np.random.multivariate_normal(initmeans[K-1][component[i],:], initcovs[K-1][component[i],:,:])

err_x, err_y, err_z, err_u, err_v, err_w = phaseP
#uvwcov_gal = np.zeros_like(uvw_gal) #use this when not using astrometry coordinate transformation


# In[ ]:

# Transform the simulated positions and velocities into astrometric observables and radial velocity in ICRS
print "Transformation to astrometric observables.."
x_icrs, y_icrs, z_icrs = galtoicrs.transformCartesianCoordinates(x_gal, y_gal, z_gal)
u_icrs, v_icrs, w_icrs = galtoicrs.transformCartesianCoordinates(uvw_gal[:,0], uvw_gal[:,1], uvw_gal[:,2])
alpha, delta, parallax, mura, mudec, vrad = phaseSpaceToAstrometry(x_icrs, y_icrs, z_icrs, u_icrs, v_icrs, w_icrs) #rad, rad, mas, mas/yr, mas/yr, km/s
astr_true = np.column_stack((alpha, delta, parallax, mura, mudec, vrad))

plt.hist(alpha)
plt.show()
plt.hist(delta)
plt.show()
plt.hist(parallax)
plt.show()
plt.hist(mura)
plt.show()
plt.hist(mudec)
plt.show()

print "Simulate measured values.."
#Simulate measured values using measurement errors
measurement_error = np.array([0.,0.,0.3,1.,1.,0])
measurement_covar = np.diag(measurement_error**2) #measurement errors
astr_measured = np.empty([N,6])
for i in range(N):
    astr_measured[i] = np.random.multivariate_normal(astr_true[i,:], measurement_covar, 1)
errors_measured = np.tile(measurement_error, (N,1))
corr_measured = np.tile(np.tile([0.],15), (N,1))

arr_astr = np.column_stack((astr_measured, errors_measured, corr_measured))
labels = ('ra', 'dec', 'parallax', 'pmra', 'pmdec', 'vrad', 'ra_error', 'dec_error', 'parallax_error', 'pmra_error', 'pmdec_error', 'vrad_error', 'ra_dec_corr', 'ra_parallax_corr', 'ra_pmra_corr', 'ra_pmdec_corr', 'ra_vrad_corr', 'dec_parallax_corr', 'dec_pmra_corr', 'dec_pmdec_corr', 'dec_vrad_corr', 'parallax_pmra_corr', 'parallax_pmdec_corr', 'parallax_vrad_corr', 'pmra_pmdec_corr', 'pmra_vrad_corr', 'pmdec_vrad_corr')
t_astr=Table(arr_astr, names=labels)


# In[39]:

print "Projection.."
# Calculate the projection matrix analogous to equation (1) in Bovy et al 2009 
# (https://ui.adsabs.harvard.edu/#abs/2009ApJ...700.1794B/abstract). Note the different ordering of the
# velocity components.
t_astr_vrad0 = cp.copy(t_astr)
for i in ['ra', 'dec']:
    gfc.add_rad(t_astr_vrad0, i, u.rad)
    gfc.add_rad(t_astr, i, u.rad)

# Do not give radial velocity data to transformation
matrix.transformation(t_astr_vrad0)#, vrad_col = 'vrad')
warr_vrad0 = gfc.XD_arr(t_astr_vrad0, "w") #icrs
wcovar_vrad0 = gfc.XD_arr(t_astr_vrad0, "S") ; wcovar_vrad0[:,2,2] = 1e15
varr_vrad0 = gfc.XD_arr(t_astr_vrad0, "UVW")
proj_vrad0 = gfc.XD_arr(t_astr_vrad0, "R")

# Give radial velocity data to transformation
matrix.transformation(t_astr, vrad_col = 'vrad')
warr = gfc.XD_arr(t_astr, "w") #icrs
wcovar = gfc.XD_arr(t_astr, "S")
varr = gfc.XD_arr(t_astr, "UVW")
proj = gfc.XD_arr(t_astr, "R")


# In[40]:

fig, ax = plt.subplots(1,3, sharey=True, figsize=(9,3))
velocities = ('U', 'V', 'W') ; unit = ' (km/s)'
for v in range(len(velocities)):
    ax[v].hist(varr_vrad0[:,v], bins='auto', normed=True, facecolor='grey', histtype='stepfilled', alpha=0.3, label='varr')
    ax[v].hist(uvw_gal[:,v], bins='auto', normed=True, facecolor='black', histtype='stepfilled', alpha=0.5, label='uvwgal')
    ax[v].set_xlabel(velocities[v] + unit)
    ax[v].legend(loc='upper right', prop={'size': 6})
plt.show()


# In[41]:

# Perform XD
print "Perform XD.."
# The input to XD are the values of v_alpha*, v_delta, vrad. The other input required is the projection matrix.
# The values of v_alpha*, v_delta, vrad are obtained from (alpha, delta, parallax, mura,  mudec, vrad).
wrange = np.arange(wmin, wmax + wstep, wstep)**2.
Krange = range(Kmin, Kmax + Kstep, Kstep) #Krange from Kmin to Kmax

logL, AIC, MDL, amps_test, means_test, covs_test, bestK, bestw = gfc.perform_XD(warr, wcovar, proj, initamps, initmeans, initcovs, wrange, Krange, N)
logL_vrad0, AIC_vrad0, MDL_vrad0, amps_test_vrad0, means_test_vrad0, covs_test_vrad0, bestK_vrad0, bestw_vrad0 = gfc.perform_XD(warr_vrad0, wcovar_vrad0, proj, initamps, initmeans, initcovs, wrange, Krange, N)


# In[42]:

print "XD:"
print "logLikelihood: best K = {0}, best w = {1}".format(bestK[0], bestw[0])
print "amps: {0}, means: {1}, covs: {2}".format(amps_test[0], means_test[0], covs_test[0])
print "AIC: best K = {0}, best w = {1}".format(bestK[1], bestw[1])
print "amps: {0}, means: {1}, covs: {2}".format(amps_test[1], means_test[1], covs_test[1])
print "MDL: best K = {0}, best w = {1}".format(bestK[2], bestw[2])
print "amps: {0}, means: {1}, covs: {2}".format(amps_test[2], means_test[2], covs_test[2])


# In[43]:

def plot_normal_PDF(ax, v, x, amps, means, covs, c, l):
    pdf = np.zeros(len(x))
    for n in range(len(amps)):
        pdf += amps[n]*sp.stats.norm.pdf(x, loc=means[n,v], scale=np.sqrt(covs[n,v,v]))
    ax.plot(x, pdf, label=l, color=c, lw=.8)

def plot_hist_uvw(inita, initm, initc, a_test, m_test, c_test, uvw_data, a_test_vrad0=None, m_test_vrad0=None, c_test_vrad0=None, uvw_data_vrad0=None, *args):
    saveto = '/disks/strw9/vanweenen/mrp1/results/'
    velocities = ('U', 'V', 'W') ; unit = ' (km/s)'
    colors_vrad0 = ('red','dodgerblue','green') ; colors = ('lightcoral', 'skyblue', 'greenyellow')
    test = ('logL', 'AIC', 'MDL') ; test_vrad0 = ('logL $v_r = 0$', 'AIC $v_r = 0$', 'MDL $v_r = 0$')
        
    fig, ax = plt.subplots(len(test), len(velocities), sharey=True, figsize=(9,9), tight_layout=True)
    for v in range(len(velocities)):
        vrange = np.linspace(np.amin(uvw_data[:,v]), np.amax(uvw_data[:,v]), len(uvw_data[:,v]))
        for t in range(len(test)):
            plot_normal_PDF(ax[t,v], v, vrange, a_test[t], m_test[t], c_test[t], colors[t], test[t]) #xd fit
            if a_test_vrad0 is not None:
                plot_normal_PDF(ax[t,v], v, vrange, a_test_vrad0[t], m_test_vrad0[t], c_test_vrad0[t], colors_vrad0[t], test_vrad0[t]) #xd fit without vrad info
            plot_normal_PDF(ax[t,v], v, vrange, inita, initm, initc, 'black', 'initial') #curve of initial velocity distribution
            ax[t,v].hist(uvw_data[:,v], bins='auto', normed=True, facecolor='black', histtype='stepfilled', alpha=0.3, label='data')
            if uvw_data_vrad0 is not None:
                ax[t,v].hist(uvw_data_vrad0[:,v], bins='auto', normed=True, facecolor='grey', histtype='stepfilled', alpha=0.15, label='data $v_r = 0$')
            ax[t,v].set_xlabel(velocities[v] + unit)
            ax[t,v].legend(loc='upper right', prop={'size': 6})
    
    suptitle = 'Histogram of velocity in Cartesian coordinates'
    filename = '/hist_velocity'
    suptitle, filename = gfc.gplot.title_filename(suptitle, filename, *args)
    plt.suptitle(suptitle, y=1., fontsize=12)
    plt.savefig(saveto + filename)
    plt.show()


# In[44]:

plot_hist_uvw(initamps[K-1], initmeans[K-1], initcovs[K-1], amps_test, means_test, covs_test, varr, amps_test_vrad0, means_test_vrad0, covs_test_vrad0, varr_vrad0, K, N, meanstep, sigma, p, q)


# In[45]:

gfc.gplot.plot_XD_w_K(logL, AIC, MDL, bestK, bestw, False, Kmin, Kmax, wmin, wmax, K, N, meanstep, sigma, p, q)
gfc.gplot.plot_XD_w_K(logL_vrad0, AIC_vrad0, MDL_vrad0, bestK_vrad0, bestw_vrad0, True, Kmin, Kmax, wmin, wmax, K, N, meanstep, sigma, p, q)


# In[ ]:

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

def contourplot(array, ax, bins):
    for i in range(2):
        for j in range(2):
            if j != i + 1:
                H, xedges, yedges = np.histogram2d(array[:,j], array[:,i+1], bins)
                cs = ax[i+j].contourf(xedges[:-1], yedges[:-1], H.T, bins, cmap = 'binary')

def set_axes_3velocities(ax):  
    velocities = ('U', 'V', 'W') ; unit = ' (km/s)'
    xpos = ('top', 'bottom') ; ypos = ('left', 'right')
    for i in range(2):
        for j in range(2):
            if j != i + 1:
                ax[i+j].set_xlim(-130,120) ; ax[i+j].set_ylim(-70,70)
                ax[i+j].set_xlabel(velocities[j] + unit) ; ax[i+j].set_ylabel(velocities[i+1] + unit)
                ax[i+j].xaxis.set_label_position(xpos[i]) ; ax[i+j].yaxis.set_label_position(ypos[j])
    ax[0].xaxis.tick_top()
    ax[0].set_xticks((-100, -50, 0, 50, 100))
    ax[0].set_yticks((-100, -50, 0, 50))
    ax[0].set_ylim(-120,60)
    ax[2].yaxis.tick_right()

def plot_Gaussian_comps(amps, means, covs, t_uvw, total=False, *args):
    saveto = '/disks/strw9/vanweenen/mrp1/results/'
    filename = '/PDF'
    if total:
        filename += '_total'
        title = "Velocity distribution"
        ext = ((-140,140),(-130,130),(-72,72))
        PDFs = map(gfc.pdf.multivariate, means, covs, amps)
        evalxyz = gfc.pdf.eval_total_PDF(PDFs, [ext[0], ext[1], ext[2]])
        evaltot = (evalxyz.sum(2), evalxyz.sum(1), evalxyz.sum(0))
        levels = np.logspace(-6.,-2.7,10)
    else:
        line = ("xy", "xz", "yz") 
        title = "Gaussian components of velocity distribution"
    title, filename = gfc.gplot.title_filename(title, filename, *args)
    fig = plt.figure(figsize=(8,8))
    plt.suptitle(title)
    ax = [plt.subplot(2,1,1), plt.subplot(2,2,3), plt.subplot(2,2,4)]
    if total:
        for i in range(2):
            for j in range(2):
                if j != i +1:
                    ax[i+j].contour(evaltot[i+j], extent = [ext[j][0], ext[j][1], ext[i+1][0], ext[i+1][1]], cmap = plt.cm.viridis, levels=levels)
    else:
        for i in range (3):
            for a, m, c in zip(amps, means, covs):
                gfc.gplot.draw_PDF_ellipse(ax[i], a, m, c, line[i], edgecolor="0.4")
    contourplot(t_uvw, ax, 10)
    set_axes_3velocities(ax)    
    KnownGroups(ax[0])
    fig.savefig(saveto + filename)
    plt.show()

plot_Gaussian_comps(amps_test[0], means_test[0], covs_test[0], uvwarr_icrs, False)
plot_Gaussian_comps(amps_test[0], means_test[0], covs_test[0], uvwarr_icrs, True)


# In[ ]:

pdfs_XD = map(gfc.pdf.multivariate, means_test[0], covs_test[0], amps_test[0])
for i,row in enumerate(t_astr_vrad0):
    vrad_predicted = gfc.radial_velocity_distribution(pdfs_XD, row["ra_rad"], row["dec_rad"], row["parallax"], row["pmra"], row["pmdec"], row["C"], t_astr['vrad_error'])
    print vrad_predicted


# In[ ]:




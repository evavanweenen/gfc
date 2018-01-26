from extreme_deconvolution import extreme_deconvolution as xd
from gfc import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import Angle
import pygaia.astrometry.vectorastrometry as pg
import scipy as sp

from matplotlib import rc
rc('font',**{'family':'serif', 'serif':['Times']})
rc('text', usetex=True)

parser = ArgumentParser()
parser.add_argument("save_folder", help = "Folder in which results will be saved")
parser.add_argument("--K", help = "Amount of components of normal distribution", default = 1, type= int)
parser.add_argument("--Rmax", help = "Maximum distance simulated stars", default=500, type=float)
parser.add_argument("--L", help = "How many stars to simulate", default=1e4, type = int)
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()
args.L = int(args.L)

if args.verbose:
    print "Finished parsing arguments"

Rmin = 1. #pc
print "K=", args.K ; print "L=", args.L

def spherical_to_carthesian(r, theta, phi):
    """
    Transform spherical to Carthesian coordinates
    Parameters
        r (distance units)
        theta (rad)
        phi (rad)
    Return
        x (distance units)
    """
    x = np.empty([args.L, 3])
    x[:,0] = r * np.sin(theta) * np.cos(phi)
    x[:,1] = r * np.sin(theta) * np.sin(phi)
    x[:,2] = r * np.cos(theta)
    return x

def simulate_x():
    """
    Generate homogeneous distribution in spherical coordinates for x
    Parameters
        phi (rad)
        theta (rad)
        r (pc)
    Return
        x (pc)
    """
    np.random.seed(1)
    phi = np.random.uniform(0, 2*np.pi, args.L)
    theta = np.arccos(np.random.uniform(-1, 1, args.L))
    r = np.random.uniform(Rmin**3, args.Rmax**3, args.L)**(1./3.) #pc
    x = spherical_to_carthesian(r, theta, phi)
    return x

def simulate_chosen_amps():
    """
    Simulate chosen amps
    """
    amps = np.tile([1./args.K], args.K)
    print "amps", amps
    amps_xd = amps * 100.
    np.save(args.save_folder + '/initial_values/initial_amps_K{K}.npy'.format(K=args.K), amps_xd)
    return amps

def simulate_random_amps():
    """
    Generate pseudo-random amplitudes (the number of stars in each component K)
    Return
        amps_prob (no units)
    """
    np.random.seed(0)
    amps = np.empty([args.K]) #these are the real amplitudes and should be used with xd
    rest = 1
    for i in range(args.K-1):
        amps[i] = np.random.uniform(0,rest)
        rest = rest - amps[i]
    amps[args.K-1] = rest
    print "amps", amps
    amps_xd = amps * 100.
    print "amps_xd", amps_xd
    np.save(args.save_folder + '/initial_values/initial_amps_K{K}.npy'.format(K=args.K), amps_xd)
    return amps

def simulate_chosen_means():
    """
    Generate chosen means
    Based on fact that average velocity of star is 100 km/s
    Return
        means (km/s)
    """
    means = np.empty([args.K, 3])
    for i in range(args.K):
        means[i,:] = [10.* i, 10.*i, 10.*i]
    print "means", means
    np.save(args.save_folder + '/initial_values/initial_means_K{K}.npy'.format(K=args.K), means)
    return means
    

def simulate_random_means():
    """
    Generate pseudo-random means, different for each component
    Based on fact that average velocity of star is 100 km/s
    Return
        means (km/s)
    """
    np.random.seed(2)
    means = np.empty([args.K, 3])
    means[:,0] = np.random.uniform(0.,200., args.K)
    means[:,1] = np.random.uniform(0.,200., args.K)
    means[:,2] = np.random.uniform(0.,200., args.K)
    print "means", means
    np.save(args.save_folder + '/initial_values/initial_means_K{K}.npy'.format(K=args.K), means)
    return means

def simulate_chosen_sigma():
    sigma2 = np.tile([.1,.1,.1], (args.K, 1))
    print "sigma2", sigma2
    return sigma2

def simulate_random_constant_sigma2():
    """
    Generate pseudo-random covariances that do not vary per component
    Return
        sigma2 (km/s)^2
    """
    np.random.seed(3)
    sigma2_i = np.random.uniform(0.,40.,3)
    sigma2 = np.tile(sigma2_i, (args.K,1))
    print "sigma2", sigma2
    return sigma2

def simulate_random_varying_sigma2():
    """
    Generate pseudo-random covariances, different for each component
    Return
        sigma2 (km/s)^2
    """
    np.random.seed(3)
    sigma2 = np.empty([args.K, 3])
    sigma2[:,0] = np.random.uniform(0.,40., args.K)
    sigma2[:,1] = np.random.uniform(0.,40., args.K)
    sigma2[:,2] = np.random.uniform(0.,40., args.K)
    print "sigma2", sigma2
    return sigma2
    
def sigma_to_cov(sigma2):
    """
    Create covariance matrix of velocity distribution with sigma2
    Parameters
        sigma2 (km/s)^2
    Return
        covs (km/s)^2
    """
    covs = np.zeros((args.K,3,3))
    covs[:,0,0] = sigma2[:,0]
    covs[:,1,1] = sigma2[:,1]
    covs[:,2,2] = sigma2[:,2]
    print "covs", covs
    np.save(args.save_folder + '/initial_values/initial_covs_K{K}.npy'.format(K=args.K), covs)
    return covs
      

def simulate_v(random, vary_sigma2 = False):
    """
    Generate normal distribution in Carthesian coordinates for velocity consisting of K components
    Parameters
        means (km/s)
        covs (km/s)^2
    Return
        v (km/s)
    """
    print "Random = ", random
    if random == False:
        amps = simulate_chosen_amps()
        means = simulate_chosen_means() #array of K by 3 means
        covs = sigma_to_cov(simulate_chosen_sigma())
    else:
        amps = simulate_random_amps()
        means = simulate_random_means()
        if vary_sigma2:
            covs = sigma_to_cov(simulate_random_varying_sigma2())
        else:
            covs = sigma_to_cov(simulate_random_constant_sigma2())
    v = np.empty([args.L,3])
    component = np.random.choice(args.K, args.L, p=amps) #choose component
    for i in range(args.L):
        v[i] = np.random.multivariate_normal(means[component[i],:], covs[component[i],:,:])
    return v, means, covs, amps

def axisnames(array, i, j, ax):
    if np.array_equal(array, x):
        names = ['X', 'Y', 'Z']
        unit = ' (pc)'
    elif np.array_equal(array, v):
        names = ['U', 'V', 'W']
        unit = ' (km/s)'
    ax.set_xlabel(names[i] + unit) 
    ax.set_ylabel(names[j] + unit)

def contourplot(array, i, j, ax, bins):
    H, xedges, yedges = np.histogram2d(array[:,i], array[:,j], bins)
    cs = ax.contourf(xedges[:-1], yedges[:-1], H, bins, cmap = plt.cm.binary)
    axisnames(array, i, j, ax)

def scatterplot(array, i, j, ax):
    ax.scatter(array[:,i], array[:,j], s=.3, c='black', linewidth = 0)
    axisnames(array, i, j, ax)

def plot(array, subj, plottype = 'contour', bins = 15):
    ax0 = plt.subplot(2,1,1)
    ax1 = plt.subplot(2,2,3)
    ax2 = plt.subplot(2,2,4)
    if plottype == 'contour':
        contourplot(array, 0, 1, ax0, bins)
        contourplot(array, 0, 2, ax1, bins)
        contourplot(array, 1, 2, ax2, bins)
    elif plottype == 'scatter':
        scatterplot(array, 0, 1, ax0)
        scatterplot(array, 0, 2, ax1)
        scatterplot(array, 1, 2, ax2)
    plt.tight_layout()
    plt.savefig(args.save_folder + '/simulated_data_{subj}_K{K}_L{L}.png'.format(subj = subj, K = args.K, L = args.L))
    #plt.show()
    
def hist_v(i, v, mean, cov, amp):
    names = ['U', 'V', 'W']
    unit = ' (km/s)'
    x = np.linspace(np.amin(v[:,i]), np.amax(v[:,i]), args.L)
    plt.figure()
    plt.title('Histogram of velocity in Carthesian coordinates')
    for k in range(args.K):
        plt.plot(x, amp[k]*sp.stats.norm.pdf(x, loc=mean[k,i], scale=np.sqrt(cov[k,i,i])), color='r', lw=1)
    plt.hist(v[:,i], facecolor='green', normed=True, histtype='stepfilled', alpha=0.3)
    plt.xlabel(names[i] + unit)
    plt.savefig(args.save_folder + '/hist_K{K}_velocity_'.format(K = args.K) + names[i] + '.png')
    #plt.show()
        
def hist_astr_err(i, distr):
    names = ['Right Ascension', 'Declination', 'Parallax', 'Proper motion of right ascension', 'Proper motion of declination']
    units = [' (mas)', ' (mas)', ' (mas)', ' (mas/yr)', ' (mas/yr)']
    plt.figure()
    plt.title('Histogram of measurement error in astrometric coordinates')
    x = np.linspace(np.amin(distr[:,i]), np.amax(distr[:,i]), args.L)
    plt.plot(x, sp.stats.norm.pdf(x, loc=0, scale=1), color='red', lw=1)
    plt.hist(distr[:,i], facecolor='green', normed=True, histtype='stepfilled', alpha=0.3)
    plt.xlabel(names[i] + units[i])
    plt.savefig(args.save_folder + '/hist_K{K}_astr_err_'.format(K = args.K) + names[i] + '.png')
    #plt.show()

"""
Simulate homogeneous distribution for x in spherical coordinates
Simulate normal distribution for v consisting of K components
"""
#Carthesian coordinates: x (pc), v (km/s)
x = simulate_x() #in pc
v, means_v, covs_v, amps_v = simulate_v(random = False) #in km/s

"""
Plot simulated data
"""
#Visualisation of data
plot(x, 'x')
plot(v, 'v')

r = np.sqrt(x[:,0]**2+x[:,1]**2+x[:,2]**2)
#a = np.linspace(Rmin, args.Rmax, args.L)
plt.figure()
#plt.plot(a, (a)**(1./3.), color='r', lw=1)
#plt.plot(a, ((3.)**(a/((args.Rmax-Rmin-499./4.)))-1)/(args.Rmax - Rmin), color='r', lw=1)
#plt.figure()
plt.hist(r, facecolor='green', normed=True, histtype='stepfilled', alpha = 0.3)
plt.xlabel(r"Radius (pc)")
plt.savefig(args.save_folder+'/hist_r.png')
#plt.show()

#Check for normal distributions
hist_v(0, v, means_v, covs_v, amps_v)
hist_v(1, v, means_v, covs_v, amps_v)
hist_v(2, v, means_v, covs_v, amps_v)

#measurement_error = (astr_measured - astr_true)/errors_tot
#hist_astr_err(2, measurement_error)
#hist_astr_err(3, measurement_error)
#hist_astr_err(4, measurement_error)

"""
Write x,v data to file
"""
arr_phasespace = np.column_stack((x,v))
labels_phasespace = ('x1','x2','x3','w1','w2','w3')
t_phasespace = Table(arr_phasespace, names=labels_phasespace)
print "t", t_phasespace
io.write_table_without_arrays(t_phasespace, args.save_folder+'/simulated_K{K}/simulated_data_xv_K{K}_L{L}.npy'.format(K = args.K, L = args.L))



"""
Transform to Astrometry coordinates
"""
#phi (rad->mas), theta (rad->mas), parallax (mas), muphistar (mas/yr), mutheta (mas/yr), vrad (km/s)
phi, theta, parallax, muphistar, mutheta, vrad = pg.phaseSpaceToAstrometry(x[:,0],x[:,1],x[:,2],v[:,0],v[:,1],v[:,2])
phi = Angle(phi, u.radian)
theta = Angle(theta, u.radian)
astr_true = np.column_stack((phi.mas, theta.mas, parallax, muphistar, mutheta)) #leave v_rad out, true values

#Measurement errors and covariances: phi_error (mas), theta_error (mas), parallax_error(mas), muphistar_error (mas/yr), mutheta_error (mas/yr)
#Take 0.3 mas for parallax and 1 mas/yr for mustar and mutheta as true variances
errors_astr = np.array([0.,0.,0.3,1.,1.])
covs_astr = np.diag(errors_astr**2)

#Simulate measured values of astrometry
astr_measured = np.empty([args.L,5])
for i in range(args.L):
    astr_measured[i] = np.random.multivariate_normal(astr_true[i,:], covs_astr, 1)

covs_tot_astr = np.tile(np.tile([0.],10), (args.L, 1)) #nondiagonal values of array
errors_tot_astr = np.tile(errors_astr, (args.L,1))

"""
Write data to table and file
"""
arr_astr = np.column_stack((astr_measured, errors_tot_astr, covs_tot_astr))
labels = ('ra', 'dec', 'ra_error', 'dec_error', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'ra_dec_corr', 'ra_parallax_corr', 'ra_pmra_corr', 'ra_pmdec_corr', 'dec_parallax_corr', 'dec_pmra_corr', 'dec_pmdec_corr', 'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr')

t_astr=Table(arr_astr, names=labels)
print "t_astr", t_astr

io.write_table_without_arrays(t_astr, args.save_folder+'/simulated_K{K}/simulated_data_astr_K{K}_L{L}.npy'.format(K = args.K, L = args.L))

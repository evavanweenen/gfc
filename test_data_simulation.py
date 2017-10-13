from gfc import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.table import Table
from pygaia.astrometry.vectorastrometry import phaseSpaceToAstrometry

from matplotlib import rc
rc('font',**{'family':'serif', 'serif':['Times']})
rc('text', usetex=True)

parser = ArgumentParser()
parser.add_argument("save_folder", help = "Folder in which results will be saved")
parser.add_argument("--K", help = "Amount of components of normal distribution", default = 1, type= int)
parser.add_argument("--errx", help = "Error in x", default=[.01,.01,.01], type=list)
parser.add_argument("--errv", help = "Error in v", default=[2,4,6], type=list)
parser.add_argument("--L", help = "How many stars to simulate", default=1e4, type = int)
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()
args.L = int(args.L)

if args.verbose:
    print "Finished parsing arguments"

Rmin = 1. #pc
Rmax = 500. #pc

print "K=", args.K ; print "L=", args.L

def spherical_to_carthesian(r, theta, phi):
    """
    Transform spherical to Carthesian coordinates
    """
    x = np.empty([args.L, 3])
    x[:,0] = r * np.sin(theta) * np.cos(phi)
    x[:,1] = r * np.sin(theta) * np.sin(phi)
    x[:,2] = r * np.cos(theta)
    return x

def simulate_x():
    """
    Generate homogeneous distribution in spherical coordinates for x
    """
    np.random.seed(1)
    phi = np.random.uniform(0, 2*np.pi, args.L)
    theta = np.arccos(np.random.uniform(-1, 1, args.L))
    r = np.random.uniform(Rmin**3, Rmax**3, args.L)**(1./3.)
    x = spherical_to_carthesian(r, theta, phi)
    return x

def simulate_amps():
    """
    Generate pseudo-random amplitudes (the number of stars in each component K)
    """
    np.random.seed(0)
    rest = 1.
    amps = []
    for i in range(args.K-1):
        amp_i = np.random.uniform(0,rest)
        rest = rest - amp_i
        amps.append(amp_i)
    amps.append(rest)
    print "amps", amps
    np.save(args.save_folder + '/simulated_K{K}/initial_amps.npy'.format(K=args.K), amps)
    return amps

def simulate_constant_means():
    """
    Generate pseudo-random means that do not vary per component
    """
    np.random.seed(2)
    mean_i = np.random.uniform(0.,100.,3)
    means = np.tile(mean_i, (args.K,1))
    print "means", means
    np.save(args.save_folder + '/simulated_K{K}/initial_means.npy'.format(K=args.K), means)
    return means
    

def simulate_varying_means():
    """
    Generate pseudo-random means, different for each component
    """
    np.random.seed(2)
    means = np.empty([args.K, 3])
    means[:,0] = np.random.uniform(0.,100., args.K)
    print means[:,0]
    means[:,1] = np.random.uniform(0.,100., args.K)
    means[:,2] = np.random.uniform(0.,100., args.K)
    print "means", means
    np.save(args.save_folder + '/simulated_K{K}/initial_means.npy'.format(K=args.K), means)
    return means

def simulate_constant_sigma2():
    """
    Generate pseudo-random covariances that do not vary per component
    """
    np.random.seed(3)
    sigma2_i = np.random.uniform(0.,40.,3)
    sigma2 = np.tile(sigma2_i, (args.K,1))
    print "sigma2", sigma2
    return sigma2

def simulate_varying_sigma2():
    """
    Generate pseudo-random covariances, different for each component
    """
    np.random.seed(3)
    sigma2 = np.empty([args.K, 3])
    sigma2[:,0] = np.random.uniform(0.,40., args.K)
    sigma2[:,1] = np.random.uniform(0.,40., args.K)
    sigma2[:,2] = np.random.uniform(0.,40., args.K)
    print "sigma2", sigma2
    return sigma2
    
def sigma_to_cov(sigma2):
    covs = np.zeros((args.K,3,3))
    covs[:,0,0] = sigma2[:,0]
    covs[:,1,1] = sigma2[:,1]
    covs[:,2,2] = sigma2[:,2]
    print "covs", covs
    np.save(args.save_folder + '/simulated_K{K}/initial_covs.npy'.format(K=args.K), covs)
    return covs
      

def simulate_v(vary_means = True, vary_sigma2 = False):
    """
    Generate normal distribution in Carthesian coordinates for velocity consisting of K components
    """
    amps = simulate_amps()
    if vary_means:
        means = simulate_varying_means() #array of K by 3 means
    else:
        means = simulate_constant_means()
    if vary_sigma2:
        covs = sigma_to_cov(simulate_varying_sigma2()) #array of K by 3 by 3 covs
    else:
        covs = sigma_to_cov(simulate_constant_sigma2())
    v = np.array([]).reshape(0,3)
    rest = args.L
    for k in range(args.K):
        l = int(amps[k] * args.L)
        if k == (args.K - 1): #if something goes wrong with rounding off l
            l = rest
        rest -= l
        v_k = np.random.multivariate_normal(means[k,:], covs[k,:,:], l)
        v = np.append(v,v_k, axis=0)
        print "K:", k+1  
    return v


"""
Simulate homogeneous distribution for x in spherical coordinates and transform to Carthesian coordinates
Simulate normal distribution for v consisting of K components in Carthesian coordinates
"""
x = simulate_x()
v = simulate_v()

#print r[r<1.]
#print x[np.sqrt(x[:,0]**2+x[:,1]**2+x[:,2]**2)<1.]
print v

"""
Plot simulated data
"""
fig = plt.figure(0)
fig.suptitle("Simulated data with K = {K} and L = {L}".format(K = args.K, L = args.L))
ax1 = fig.add_subplot(1,2,1, projection='3d') ; ax2 = fig.add_subplot(1,2,2, projection='3d')
ax1.scatter(x[:,0], x[:,1], x[:,2], s=.3, c='black', linewidth = 0)
ax2.scatter(v[:,0], v[:,1], v[:,2], s=.3, c='black', linewidth = 0)
ax1.set_xlabel("x") ; ax2.set_xlabel("U")
ax1.set_ylabel("y") ; ax2.set_ylabel("V")
ax1.set_zlabel("z") ; ax2.set_zlabel("W")
plt.show()
fig.savefig(args.save_folder + '/simulated_data_K{K}_L{L}.png'.format(K = args.K, L = args.L))

"""
Transform simulated data in Carthesian (phasespace) coordinates to observational (astrometry) coordinates
"""
#positions
phi, theta, parallax, muphistar, mutheta, vrad = phaseSpaceToAstrometry(x[:,0],x[:,1],x[:,2],v[:,0],v[:,1],v[:,2])
positions = np.column_stack((phi, theta, parallax, muphistar, mutheta))
print positions

#errors
e = .1 #simple error
errors = np.tile(np.tile([e],5), (args.L,1))

#covariances
c = 0. #simple covariances
covariances = np.tile(np.tile([c],10), (args.L, 1))

"""
Write data to table and file
"""
arr = np.column_stack((positions, errors, covariances))
labels = ('ra', 'dec', 'ra_error', 'dec_error', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'ra_dec_corr', 'ra_parallax_corr', 'ra_pmra_corr', 'ra_pmdec_corr', 'dec_parallax_corr', 'dec_pmra_corr', 'dec_pmdec_corr', 'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr')
print len(arr[:,0]), len(arr[0,:])
print len(labels)

t=Table(arr, names=labels)

io.write_table_without_arrays(t, args.save_folder+'/simulated_K{K}/simulated_data_K{K}_L{L}.npy'.format(K = args.K, L = args.L))

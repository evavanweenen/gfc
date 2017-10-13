import gfc
from gfc import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from pygaia.astrometry.vectorastrometry import phaseSpaceToAstrometry

parser = ArgumentParser()
parser.add_argument("save_folder", help = "Folder in which results will be saved")
parser.add_argument("--R", help = "Maximum radius of position of simulated stars", default = 15., type=float)
parser.add_argument("--K", help = "Amount of components of normal distribution", default = 1, type= int)
parser.add_argument("--errx", help = "Error in x", default=[.01,.01,.01], type=list)
parser.add_argument("--errv", help = "Error in v", default=[2,4,6], type=list)
parser.add_argument("--L", help = "How many stars to simulate", default=1e2, type = int)
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()
args.L = int(args.L)

if args.verbose:
    print "Finished parsing arguments"

def simulate_v(cov_v):
    """
    Generate normal distribution in Carthesian coordinates for velocity consisting of K components
    Here, mean is dependent on K (to make varying easy)
    Covariances stay the same.
    """
    mean_tot = np.array([]).reshape(0,3)
    cov_tot = np.array([]).reshape(0,3,3)
    v_tot = np.array([]).reshape(0,3)
    for i in range(args.K):
        mean_v = np.array([10*i,10*i,10*i])
        length_i = args.L/args.K
        if i == (args.K - 1):
            length_i = args.L - i * length_i
        print "length", length_i
        v_i = np.random.multivariate_normal(mean_v, cov_v, length_i)
        v_tot = np.append(v_tot,v_i, axis=0)
        #mean_tot = np.append(mean_tot, [mean_v], axis=0)
        #cov_tot = np.append(cov_tot, [cov_v], axis=0)
        i+= 1
        print "K:", i
        
    return v_tot#, mean_tot, cov_tot

"""
Simulate data in Carthesian coordinates
"""
#Homogeneous distribution for x in radians
phi = np.random.uniform(0, 2*np.pi, args.L)
costheta = np.random.uniform(-1, 1, args.L)
theta = np.arccos(costheta)
r = args.R * np.random.uniform(0, 1, args.L)**(1./3.)

x=np.empty([args.L, 3])
x[:,0] = r * np.sin(theta) * np.cos(phi)
x[:,1] = r * np.sin(theta) * np.sin(phi)
x[:,2] = r * np.cos(theta)

#3D Normal distribution for v
cov_v = np.array([[(args.errv[0])**2,0,0], [0,(args.errv[1])**2,0], [0,0,(args.errv[2])**2]])  # diagonal covariance
v = simulate_v(cov_v) #x,y,z
plt.scatter(v[:,0], v[:,1])
plt.show()

#np.save(args.save_folder + 'initial_means.npy', initial_mean)
#np.save(args.save_folder + 'initial_covs.npy', initial_cov)


"""
Transform simulated data in Carthesian (phasespace) coordinates to observational (astrometry) coordinates
"""
#positions
phi, theta, parallax, muphistar, mutheta, vrad = phaseSpaceToAstrometry(x[:,0],x[:,1],x[:,2],v[:,0],v[:,1],v[:,2])
positions = np.column_stack((phi, theta, parallax, muphistar, mutheta))

#errors
phi_err, theta_err, parallax_err, muphistar_err, mutheta_err, vrad_err = phaseSpaceToAstrometry(args.errx[0], args.errx[1], args.errx[2], args.errv[0], args.errv[1], args.errv[2])
#very ugly manner
phi_err = [phi_err] * args.L
theta_err = [theta_err] * args.L
parallax_err = [parallax_err] * args.L
muphistar_err = [muphistar_err] * args.L
mutheta_err = [mutheta_err] * args.L
vrad_err = [vrad_err] * args.L
errors = np.column_stack((phi_err, theta_err, parallax_err, muphistar_err, mutheta_err))

#covariances
phitheta_cov = [0] * args.L
phiparallax_cov = [0] * args.L
phimustar_cov = [0] * args.L
phimutheta_cov = [0] * args.L
thetaparallax_cov = [0] * args.L
thetamustar_cov = [0] * args.L
thetamutheta_cov = [0] * args.L
parallaxmustar_cov = [0] * args.L
parallaxmutheta_cov = [0] * args.L
mustarmutheta_cov = [0] * args.L
covariances = np.column_stack((phitheta_cov, phiparallax_cov, phimustar_cov, phimutheta_cov, thetaparallax_cov, thetamustar_cov, thetamutheta_cov, parallaxmustar_cov, parallaxmutheta_cov, mustarmutheta_cov))


"""
Write data to table and file
"""
arr = np.column_stack((positions, errors, covariances))
labels = ('ra', 'dec', 'ra_error', 'dec_error', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'ra_dec_corr', 'ra_parallax_corr', 'ra_pmra_corr', 'ra_pmdec_corr', 'dec_parallax_corr', 'dec_pmra_corr', 'dec_pmdec_corr', 'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr')
print len(arr[:,0]), len(arr[0,:])
print len(labels)

t=Table(arr, names=labels)

gfc.io.write_table_without_arrays(t, args.save_folder+'/simulated_data_K{K}.npy'.format(K = args.K))

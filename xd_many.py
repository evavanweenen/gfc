"""
Fit a model to TGAS data using XD.
"""
import gfc

from gfc import ArgumentParser
parser = ArgumentParser()
parser.add_argument("data_file", help = "File that contains the data")
parser.add_argument("save_folder", help = "Folder in which results will be saved")
parser.add_argument("--init_amps", help = "File with initial estimates of amplitudes", default = "Bovy_parameters/Bovy_amps.npy")
parser.add_argument("--init_means", help = "File with initial estimates of means", default = "Bovy_parameters/Bovy_means.npy")
parser.add_argument("--init_covs", help = "File with initial estimates of covariances", default = "Bovy_parameters/Bovy_covs.npy")
parser.add_argument("--Kmin", help = "Lowest value of K to test", default = 4, type = int)
parser.add_argument("--Kmax", help = "Highest value of K to test", default = 20, type = int)
parser.add_argument("--Kstep", help = "Step between K values", default = 1, type = int)
parser.add_argument("--wmin", help = "Lowest value of SQRT(w) to test", default = 0.5, type = float)
parser.add_argument("--wmax", help = "Highest value of SQRT(w) to test", default = 3.5, type = float)
parser.add_argument("--wstep", help = "Step between SQRT(w) values", default = 0.25, type = float)
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

if args.verbose:
    print "Finished parsing arguments"

"""
Read initial amplitudes
"""
initial_amps = gfc.io.load(args.init_amps)
initial_means = gfc.io.load(args.init_means)
initial_covs = gfc.io.load(args.init_covs)

if args.verbose:
    print "Loaded initial estimates for Gaussian parameters"


"""
Create w and K array
"""    
wrange = gfc.np.arange(args.wmin, args.wmax + args.wstep, args.wstep)**2.
Krange = range(args.Kmin, args.Kmax + args.Kstep, args.Kstep)

"""
Create separate folder for each different model (so for each w and K combination.
Keep in mind that this is not the best solution for the w-K likelihood plot.
"""
for w in wrange:
    f_w = "{f}/w_{w}".format(f = args.save_folder, w = w)
    if not gfc.isdir(f_w):
        gfc.mkdir(f_w)
    for K in Krange:
        f = "{f}/w_{w}/K_{K}".format(f = args.save_folder, w = w, K = K)
        if not gfc.isdir(f):
            gfc.mkdir(f)

if args.verbose:
    print "Finished creating folders"


"""
Read data
"""
time_before_loading = gfc.time() #times how long programme is running
t = gfc.io.read_csv(args.data_file)
gfc.remove_unused_columns(t)
time_after_loading = gfc.time()

if args.verbose:
    print "Finished loading data in {0:.1f} seconds".format(time_after_loading - time_before_loading)


"""
Calculate vectors and covariance matrix
"""
time_before_matrices = gfc.time()
gfc.add_rad(t)
gfc.matrix.add_w(t)
gfc.matrix.add_A(t)
gfc.matrix.add_R(t)
gfc.tgas.add_C(t)
gfc.matrix.add_Q(t)
gfc.matrix.add_S(t)
time_after_matrices = gfc.time()
if args.verbose:
    print "Added w, A, R, C, Q, S in {0:.1f} seconds".format(time_after_matrices - time_before_matrices)


"""
XD cannot handle csv format, so convert to numpy tables
"""
warr = gfc.XD_arr(t, "w1", "w2", "w3")
wcov = gfc.XD_arr(t, "S") ; wcov[:, 0, 0] = 1e15 #hardcode: pretend first component (vrad) is 0, and make covariance very high so it is ignored
proj = gfc.XD_arr(t, "R")


"""
Perform XD and write to file
"""
for K in Krange:
    for w in wrange:
        f = "{f}/w_{w}/K_{K}".format(f = args.save_folder, w = w, K = K)
        print "K = {0} ; w = {1}".format(K, w)
        amps_xd, means_xd, covs_xd, L = gfc.XD(warr, wcov, initial_amps[:K], initial_means[:K], initial_covs[:K], projection = proj, w = w) #extreme deconvolution, if data is simulated and it does not work out, something goes wrong at proj
        print >> open("{0}/L".format(f), 'w'), L
        gfc.io.save_PDFs(amps_xd, means_xd, covs_xd, f)

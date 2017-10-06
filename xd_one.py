"""
Fit a model to TGAS data using XD.
Done for one model (so one k and W), and get amps, means, covs and likelihoods for one k and one W value.
"""
import gfc

from gfc import ArgumentParser
#when running programm, add options from terminal
parser = ArgumentParser()
parser.add_argument("data_file", help = "File that contains the data")
parser.add_argument("save_folder", help = "Folder in which results will be saved")
parser.add_argument("-K", help = "K", default = 15, type = int)
parser.add_argument("-w", help = "w", default = 4.0, type = float)
parser.add_argument("--init_amps", help = "File with initial estimates of amplitudes", default = "Bovy_parameters/Bovy_amps.npy") #niet gebruiken, dus standaard bovy gebruiken
parser.add_argument("--init_means", help = "File with initial estimates of means", default = "Bovy_parameters/Bovy_means.npy") #niet gebruiken 
parser.add_argument("--init_covs", help = "File with initial estimates of covariances", default = "Bovy_parameters/Bovy_covs.npy") #niet gebruiken
parser.add_argument("-v", "--verbose", action = "store_true") #always use
args = parser.parse_args()

if args.verbose:
    print "Finished parsing arguments"

#load files with initial estimates
initial_amps = gfc.io.load(args.init_amps)
initial_means = gfc.io.load(args.init_means)
initial_covs = gfc.io.load(args.init_covs)
K = args.K
w = args.w

if args.verbose:
    print "Loaded initial estimates for Gaussian parameters"

#read data (in csv) (data is for example gaia data)
time_before_loading = gfc.time() #timen hoe lang dingen duren
t = gfc.io.read_csv(args.data_file)
gfc.remove_unused_columns(t)
time_after_loading = gfc.time()

if args.verbose:
    print "Finished loading data in {0:.1f} seconds".format(time_after_loading - time_before_loading)

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

warr = gfc.XD_arr(t, "w1", "w2", "w3")
wcov = gfc.XD_arr(t, "S") ; wcov[:, 0, 0] = 1e15
proj = gfc.XD_arr(t, "R")

print "K = {0} ; w = {1}".format(K, w)
amps_xd, means_xd, covs_xd, L = gfc.XD(warr, wcov, initial_amps[:K], initial_means[:K], initial_covs[:K], projection = proj, w = w) #extreme deconvolution
print >> open("{0}/L".format(args.save_folder), 'w'), L #write loglikelihood to file
gfc.io.save_PDFs(amps_xd, means_xd, covs_xd, args.save_folder) #write 3 npy files with aps, means, covs of pdf

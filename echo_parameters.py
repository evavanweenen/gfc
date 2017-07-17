import numpy as np
import argparse
import time

parser = argparse.ArgumentParser()
parser.add_argument("xd_results_folder", help = "Folder that contains the XD results")
args = parser.parse_args()

amps_xd = np.load(args.xd_results_folder+"amplitudes.npy")
means_xd = np.load(args.xd_results_folder+"means.npy")
covs_xd = np.load(args.xd_results_folder+"covariances.npy")

print amps_xd
raw_input()
print means_xd
raw_input()
print covs_xd

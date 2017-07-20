import gfc
from gfc import ArgumentParser

parser = ArgumentParser()
parser.add_argument("xd_results_folder", help = "Folder that contains the XD results")
args = parser.parse_args()

amps_xd, means_xd, covs_xd = gfc.io.load_PDFs(args.xd_results_folder)

print amps_xd
raw_input()
print means_xd
raw_input()
print covs_xd

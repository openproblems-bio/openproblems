## VIASH START
# The code between the the comments above and below gets stripped away before 
# execution. Here you can put anything that helps the prototyping of your script.
par = {
    "input_rna": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz",
    "input_adt": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DADT%5Fumi%2Ecsv%2Egz",
    "output": "output.h5ad",
    "test": False
}
resources_dir = "../../resources/utils.py"
## VIASH END

print("Importing libraries")
import scprep

# adding resources dir to system path
# the resources dir contains all files listed in the '.functionality.resources' part of the
# viash config, amongst which is the 'utils.py' file we need.
import sys
sys.path.append(resources_dir)

# importing helper functions from common utils.py file in resources dir
from utils import create_joint_adata
from utils import filter_joint_data_empty_cells
from utils import subset_joint_data

print("Downloading expression datasets from GEO (this might take a while)") 
sys.stdout.flush()

# par["input_rna"] can be the path to a local file, or a url
rna_data = scprep.io.load_csv(
    par["input_rna"], cell_axis="col", compression="gzip", sparse=True, chunksize=1000
)
adt_data = scprep.io.load_csv(
    par["input_adt"], cell_axis="col", compression="gzip", sparse=True, chunksize=1000
)

print("Transforming into adata")
adata = create_joint_adata(rna_data, adt_data)
adata = filter_joint_data_empty_cells(adata)

adata.uns["dataset_name"] = "citeseq_cbmc"

if par["test"]:
    print("Subsetting dataset")
    adata = subset_joint_data(adata)
    adata.uns["dataset_name"] = "citeseq_cbmc_test"

print("Writing adata to file")
adata.write(par["output"], compression = "gzip")
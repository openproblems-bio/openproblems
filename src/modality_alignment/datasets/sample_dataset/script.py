print("Importing libraries")
import scprep
import pandas as pd
import numpy as np
import scipy.sparse

# adding resources dir to system path
# the resources dir contains all files listed in the '.functionality.resources' part of the
# viash config, amongst which is the 'utils.py' file we need.
import sys
sys.path.append(resources_dir)

# importing helper functions from common utils.py file in resources dir
from utils import create_joint_adata
from utils import filter_joint_data_empty_cells
from utils import subset_joint_data


rna_cells_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3271044"
    "&format=file&file=GSM3271044%5FRNA%5Fmouse%5Fkidney%5Fcell.txt.gz"
)
rna_genes_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3271044"
    "&format=file&file=GSM3271044%5FRNA%5Fmouse%5Fkidney%5Fgene.txt.gz"
)
atac_cells_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3271045"
    "&format=file&file=GSM3271045%5FATAC%5Fmouse%5Fkidney%5Fcell.txt.gz"
)
atac_genes_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3271045"
    "&format=file&file=GSM3271045%5FATAC%5Fmouse%5Fkidney%5Fpeak.txt.gz"
)

print("Downloading input files") 
sys.stdout.flush()
rna_genes = pd.read_csv(rna_genes_url, low_memory=False, index_col=0)
atac_genes = pd.read_csv(atac_genes_url, low_memory=False, index_col=1)
rna_cells = pd.read_csv(rna_cells_url, low_memory=False, index_col=0)
atac_cells = pd.read_csv(atac_cells_url, low_memory=False, index_col=0)

print("Creating joint adata object") 
keep_cells = np.intersect1d(rna_cells.index, atac_cells.index)[:200]
rna_cells = rna_cells.loc[keep_cells]
atac_cells = atac_cells.loc[keep_cells]

rna_data = scipy.sparse.csr_matrix((len(keep_cells), len(rna_genes)))
atac_data = scipy.sparse.csr_matrix((len(keep_cells), len(atac_genes)))

adata = create_joint_adata(
    rna_data,
    atac_data,
    X_index=rna_cells.index,
    X_columns=rna_genes.index,
    Y_index=atac_cells.index,
    Y_columns=atac_genes.index,
)

print("Merging obs and var") 
adata.obs = rna_cells.loc[adata.obs.index]
adata.var = rna_genes
for key in atac_cells.columns:
    adata.obs[key] = atac_cells[key]
adata.uns["mode2_varnames"] = []
for key in atac_genes.columns:
    varname = "mode2_var_{}".format(key)
    adata.uns[varname] = atac_genes[key].values
    adata.uns["mode2_varnames"].append(varname)

adata.X = scipy.sparse.csr_matrix(np.random.poisson(0.1, adata.X.shape)).astype(np.float64)
adata.obsm["mode2"] = scipy.sparse.csr_matrix(
    np.random.poisson(0.1, adata.obsm["mode2"].shape)
).astype(np.float64)

adata = filter_joint_data_empty_cells(adata)

print("Subsetting dataset")
adata = subset_joint_data(adata, n_cells = par["n_cells"], n_genes = par["n_genes"])

adata.uns["dataset_id"] = "sample_dataset_test"

print("Writing adata to file")
adata.write_h5ad(par["output"], compression = "gzip")

import anndata as ad
import pandas as pd
import numpy as np
from scipy import sparse

## VIASH START
par = {
  "input": "GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad",
  "mod1": "GEX",
  "mod2": "ATAC",
  "dataset_id": "openproblems/neurips2021_bmmc",
  "dataset_name": "BMMC (CITE-seq)",
  "dataset_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122",
  "dataset_reference": "Neurips",
  "dataset_summary": "value",
  "dataset_description": "value",
  "dataset_organism": "homo_sapiens",
  "output_mod1": "output/mod1.h5ad",
  "output_mod2": "output/mod2.h5ad"
}
meta = {
  "name": "openproblems_neurips2021_bmmc",
  "resources_dir": "/tmp/viash_inject_openproblems_neurips2021_bmmc14365472827677740971", 
}
## VIASH END

def remove_mod_col(df, mod):
  df.drop(list(df.filter(like=mod)), axis=1, inplace=True)

def remove_mod_prefix(df, mod):
  suffix = f"{mod}_"
  df.columns = df.columns.str.removeprefix(suffix)

def convert_matrix(adata):
  for key in adata:
      if isinstance(adata[key], sparse.csr_matrix):
        adata[key] = sparse.csc_matrix(adata[key])
      

print("load dataset file", flush=True)
adata = ad.read_h5ad(par["input"])

# Convert to sparse csc_matrix
convert_matrix(adata.layers)
convert_matrix(adata.obsm)

# Add is_train to obs if it is missing
if "is_train" not in adata.obs.columns:
  batch_info = adata.obs["batch"]
  batch_categories = batch_info.dtype.categories
  # From https://github.com/openproblems-bio/neurips2021_multimodal_viash/blob/75281c039ab98b459edbf52058a18597e710ed4d/src/common/datasets/process_inhouse_datasets/script.R#L14-L17
  train = ["s1d1", "s1d2", "s2d1", "s2d4", "s3d1", "s3d6", "s3d7"]
  adata.obs["is_train"] = [ "train" if x in train else "test" for x in batch_info ]

# Construct Modality datasets
print("Construct Mod datasets", flush=True)
mask_mod1 = adata.var['feature_types'] == par["mod1"]
mask_mod2 = adata.var['feature_types'] == par["mod2"]

adata_mod1 = adata[:, mask_mod1]
adata_mod2 = adata[:, mask_mod2]

# Remove other modality data from obs and var
mod1_var = pd.DataFrame(adata_mod1.var)
remove_mod_col(mod1_var, par["mod2"])
remove_mod_prefix(mod1_var, par["mod1"])
mod1_var.index.name = "feature_name"
mod1_var.reset_index("feature_name", inplace=True)
mod1_var["feature_id"] = np.where(mod1_var.gene_id.isna(), mod1_var.feature_name, mod1_var.gene_id.astype(str))
mod1_var.drop("gene_id", axis=1, inplace=True)
mod1_var.set_index("feature_id", drop=False, inplace=True)

mod1_obs = pd.DataFrame(adata_mod1.obs)
remove_mod_col(mod1_obs, par["mod2"])
remove_mod_prefix(mod1_obs, par["mod1"])

adata_mod1.var = mod1_var
adata_mod1.obs = mod1_obs

adata_mod1.uns = { key.replace(f"{par['mod1']}_", ""): value for key, value in adata.uns.items() if not key.startswith(par['mod2'])}
del adata_mod1.obsm
del adata_mod1.X

mod2_var = pd.DataFrame(adata_mod2.var)
remove_mod_col(mod2_var, par["mod1"])
remove_mod_prefix(mod2_var, par["mod2"])
mod2_var.index.name = "feature_name"
mod2_var.reset_index("feature_name", inplace=True)
mod2_var["feature_id"] = np.where(mod2_var.gene_id.isna(), mod2_var.feature_name, mod2_var.gene_id.astype(str))
mod2_var.drop("gene_id", axis=1, inplace=True)
mod2_var.set_index("feature_id", drop=False, inplace=True)

mod2_obs = pd.DataFrame(adata_mod2.obs)
remove_mod_col(mod2_obs, par["mod1"])
remove_mod_prefix(mod2_obs, par["mod2"])

adata_mod2.var = mod2_var
adata_mod2.obs = mod2_obs

adata_mod2.uns = { key.replace(f"{par['mod2']}_", ""): value for key, value in adata.uns.items() if not key.startswith(par['mod1'])}
if par["mod2"] == "ATAC":
  adata_mod2.obsm = { key.replace(f"{par['mod2']}_", ""): value for key, value in adata_mod2.obsm.items() if key.startswith(par['mod2'])}
else:
  del adata_mod2.obsm


del adata_mod2.X

print("Add metadata to uns", flush=True)
metadata_fields = [
    "dataset_id", "dataset_name", "dataset_url", "dataset_reference",
    "dataset_summary", "dataset_description", "dataset_organism"
]
for key in metadata_fields:
    if key in par:
        print(f"  Setting .uns['{key}']", flush=True)
        adata_mod1.uns[key] = par[key]
        adata_mod2.uns[key] = par[key]

print("Writing adata to file", flush=True)
adata_mod1.write_h5ad(par["output_mod1"], compression="gzip")
adata_mod2.write_h5ad(par["output_mod2"], compression="gzip")





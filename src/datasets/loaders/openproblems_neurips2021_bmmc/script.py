import anndata as ad
import pandas as pd

## VIASH START
par = {
  "input": "bmmc_multiome.decompress_gzip.h5ad",
  "mod1": "GEX",
  "mod2": "ATAC",
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
  "functionality_name": "openproblems_neurips2021_bmmc",
  "resources_dir": "/tmp/viash_inject_openproblems_neurips2021_bmmc14365472827677740971", 
}
## VIASH END

def remove_other_mod_col(df, mod):
  df.drop(list(df.filter(like=mod)), axis=1, inplace=True)

def remove_mod_prefix(df, mod):
  suffix = f"{mod}_"
  df.columns = df.columns.str.removeprefix(suffix)


print("load dataset file", flush=True)
adata = ad.read_h5ad(par["input"])

# Construct Modality datasets
print("Construct Mod datasets", flush=True)
mask_mod1 = adata.var['feature_types'] == par["mod1"]
mask_mod2 = adata.var['feature_types'] == par["mod2"]

adata_mod1 = adata[:, mask_mod1]
adata_mod2 = adata[:, mask_mod2]

# Remove other modality data from obs and var
mod1_var = pd.DataFrame(adata_mod1.var)
remove_other_mod_col(mod1_var, par["mod2"])
remove_mod_prefix(mod1_var, par["mod1"])
mod1_var.index.name = "gene_symbol"
mod1_var.reset_index("gene_symbol", inplace=True)
mod1_var.set_index("gene_id", inplace=True)

mod1_obs = pd.DataFrame(adata_mod1.obs)
remove_other_mod_col(mod1_obs, par["mod2"])
remove_mod_prefix(mod1_obs, par["mod1"])

adata_mod1.var = mod1_var
adata_mod1.obs = mod1_obs

adata_mod1.uns = { key.replace(f"{par['mod1']}_", ""): value for key, value in adata.uns.items() if not key.startswith(par['mod2'])}
del adata_mod1.obsm
del adata_mod1.X

mod2_var = pd.DataFrame(adata_mod2.var)
remove_other_mod_col(mod2_var, par["mod1"])
remove_mod_prefix(mod2_var, par["mod2"])
mod2_var.gene_id = mod2_var.index.values
mod2_var.index.name = "gene_symbol"
mod2_var.reset_index("gene_symbol", inplace=True)
mod2_var.set_index("gene_id", inplace=True)

mod2_obs = pd.DataFrame(adata_mod2.obs)
remove_other_mod_col(mod2_obs, par["mod1"])
remove_mod_prefix(mod2_obs, par["mod2"])

adata_mod2.var = mod2_var
adata_mod2.obs = mod2_obs

adata_mod2.uns = { key.replace(f"{par['mod2']}_", ""): value for key, value in adata.uns.items() if not key.startswith(par['mod1'])}
if par["mod2"] == "ATAC":
  adata_mod2.obsm = { key.replace(f"{par['mod2']}_", ""): value for key, value in adata_mod2.uns.items() if key.startswith(par['mod2'])}
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

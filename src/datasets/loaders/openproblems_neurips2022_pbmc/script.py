import anndata as ad
from scipy import sparse

## VIASH START
par = {
  "input_mod1": "cite_rna_merged.h5ad",
  "input_mod2": "cite_prot_merged.h5ad",
  "mod1": "GEX",
  "mod2": "ADT",
  "dataset_id": "openproblems/neurips2022_pbmc",
  "dataset_name": "Kaggle22 PBMC (CITE-seq)",
  "dataset_url": "https://www.kaggle.com/competitions/open-problems-multimodal/data",
  "dataset_reference": "Neurips22",
  "dataset_summary": "Neurips22 competition dataset",
  "dataset_description": "The dataset for this competition comprises single-cell multiomics data collected from mobilized peripheral CD34+ hematopoietic stem and progenitor cells (HSPCs) isolated from four healthy human donors.",
  "dataset_organism": "homo_sapiens",
  "output_mod1": "output/mod1.h5ad",
  "output_mod2": "output/mod2.h5ad"
}
meta = {
  "functionality_name": "openproblems_neurips2022_pbmc",
}
## VIASH END


def convert_matrix(adata):
  for key in adata:
      if isinstance(adata[key], sparse.csr_matrix):
        adata[key] = sparse.csc_matrix(adata[key])
      

print("load dataset modality 1 file", flush=True)
adata_mod1 = ad.read_h5ad(par["input_mod1"])

print("load dataset modality 2 file", flush=True)
adata_mod2 = ad.read_h5ad(par["input_mod2"])

# Convert to sparse csc_matrix
convert_matrix(adata_mod1.layers)
convert_matrix(adata_mod1.obsm)
convert_matrix(adata_mod2.layers)
convert_matrix(adata_mod2.obsm)


# Add is_train to obs (modality 1)
if "is_train" not in adata_mod1.obs.columns:
    split_info = adata_mod1.obs["kaggle_dataset"]
    train_sets = ["train", "test_public"]
    adata_mod1.obs["is_train"] = [ "train" if x in train_sets else "test" for x in split_info ]

# Add is_train to obs if it is missing (modality 2)
if "is_train" not in adata_mod2.obs.columns:
    split_info = adata_mod2.obs["kaggle_dataset"]
    train_sets = ["train", "test_public"]
    adata_mod2.obs["is_train"] = [ "train" if x in train_sets else "test" for x in split_info ]


# split up index in modality 1 into feature ID and feature name
adata_mod1.var['feature_id'] = [str(s).split('_')[0] for s in adata_mod1.var.index.tolist()]
adata_mod1.var['feature_name'] = [str(s).split('_')[1] for s in adata_mod1.var.index.tolist()]
adata_mod1.var.set_index('feature_id',drop=False, inplace=True)

# set feature_name (proteins have only partial ensmble IDs))
adata_mod2.var['feature_name'] = adata_mod2.var.index.tolist()
adata_mod2.var.set_index('feature_name',drop=False, inplace=True)


# remove adata.X
del adata_mod1.X
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





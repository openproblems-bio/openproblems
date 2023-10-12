import anndata as ad

## VIASH START
par = {
  "input_mod1": "resources_test/match_modalities/scicar_cell_lines/dataset_mod1.h5ad",
  "input_mod2": "resources_test/match_modalities/scicar_cell_lines/dataset_mod2.h5ad",
  "input_solution_mod1": "resources_test/match_modalities/scicar_cell_lines/solution_mod1.h5ad",
  "input_solution_mod2": "resources_test/match_modalities/scicar_cell_lines/solution_mod2.h5ad",
  "output_mod1": "output.mod1.h5ad",
  "output_mod2": "output.mod2.h5ad",
}
meta = {
    "functionality_name": "true_features"
}
## VIASH END

print("Reading input h5ad file", flush=True)
adata_mod1 = ad.read_h5ad(par["input_mod1"])
adata_mod2 = ad.read_h5ad(par["input_mod2"])

print("Storing true features", flush=True)
output_mod1 = ad.AnnData(
  obs=adata_mod1.obs[[]],
  var=adata_mod1.var[[]],
  obsm={
    "integrated": adata_mod1.obsm["X_svd"]
  },
  uns={
    "dataset_id": adata_mod1.uns["dataset_id"],
    "normalization_id": adata_mod1.uns["normalization_id"],
    "method_id": meta["functionality_name"]
  }
)
output_mod2 = ad.AnnData(
  obs=adata_mod2.obs[[]],
  var=adata_mod2.var[[]],
  obsm={
    "integrated": adata_mod2.obsm["X_svd"]
  },
  uns={
    "dataset_id": adata_mod2.uns["dataset_id"],
    "normalization_id": adata_mod2.uns["normalization_id"],
    "method_id": meta["functionality_name"]
  }
)

print("Write output to file", flush=True)
output_mod1.write_h5ad(par["output_mod1"], compression="gzip")
output_mod2.write_h5ad(par["output_mod2"], compression="gzip")

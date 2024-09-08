import anndata as ad
import harmonicalignment

## VIASH START
par = {
  "mod1" : "resources_test/common/scicar_cell_lines/dataset_mod1.h5ad",
  "mod2" : "resources_test/common/scicar_cell_lines/dataset_mod2.h5ad",
  "output" : "output.scot.h5ad",
  "n_pca_XY" : 100,
  "eigenvectors" : 100
}
meta = {
  "functionality_name" : "harmonic_alignment"
}
## VIASH END


print("Reading input h5ad file", flush=True)
adata_mod1 = ad.read_h5ad(par["input_mod1"])
adata_mod2 = ad.read_h5ad(par["input_mod2"])

print("Check parameters", flush=True)
n_eigenvectors = par["n_eigenvectors"]
n_pca_XY = par["n_pca_XY"]

if adata_mod1.layers["normalized"].shape[0] <= n_eigenvectors:
    n_eigenvectors = None
if adata_mod1.layers["normalized"].shape[0] <= n_pca_XY:
    n_pca_XY = None


print("Running Harmonic Alignment", flush=True)
ha_op = harmonicalignment.HarmonicAlignment(
    n_filters=8, n_pca_XY=n_pca_XY, n_eigenvectors=n_eigenvectors
)
ha_op.align(adata_mod1.obsm["X_svd"], adata_mod2.obsm["X_svd"])
XY_aligned = ha_op.diffusion_map(n_eigenvectors=n_eigenvectors)

print("Storing output data structures", flush=True)

adata_mod1.obsm["integrated"] = XY_aligned[: adata_mod1.obsm["X_svd"].shape[0]]
adata_mod2.obsm["integrated"] = XY_aligned[-adata_mod2.obsm["X_svd"].shape[0] :]

print("Write output to file", flush=True)
adata_mod1.uns["method_id"] = meta["functionality_name"]
adata_mod2.uns["method_id"] = meta["functionality_name"]
adata_mod1.write_h5ad(par["output_mod1"], compression = "gzip")
adata_mod2.write_h5ad(par["output_mod2"], compression = "gzip")

print("Loading dependencies")
import scanpy as sc
import harmonicalignment
import sklearn.decomposition

## VIASH START
par = {
  input = "output.h5ad",
  output = "output.scot.h5ad",
  n_svd = 100,
  n_pca_XY = 100
  eigenvectors = 100
}
resources_dir = "../../utils/"
## VIASH END

# importing helper functions from common preprocessing.py file in resources dir
import sys
sys.path.append(resources_dir)
from preprocessing import log_cpm
from preprocessing import sqrt_cpm

print("Reading input h5ad file")
adata = sc.read_h5ad(par["input"])

print("Check parameters")
n_svd = min([par["n_svd"], min(adata.X.shape) - 1, min(adata.obsm["mode2"].shape) - 1])
n_eigenvectors = par["n_eigenvectors"]
n_pca_XY = par["n_pca_XY"]

if adata.X.shape[0] <= n_eigenvectors:
    n_eigenvectors = None
if adata.X.shape[0] <= n_pca_XY:
    n_pca_XY = None

print("Normalising mode 1")
sqrt_cpm(adata)

print("Normalising mode 2")
log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")

print("Performing PCA reduction")
X_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.X)
Y_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.obsm["mode2"])

print("Running Harmonic Alignment")
ha_op = harmonicalignment.HarmonicAlignment(
    n_filters=8, n_pca_XY=n_pca_XY, n_eigenvectors=n_eigenvectors
)
ha_op.align(X_pca, Y_pca)
XY_aligned = ha_op.diffusion_map(n_eigenvectors=n_eigenvectors)

print("Storing output data structures")
adata.obsm["aligned"] = XY_aligned[: X_pca.shape[0]]
adata.obsm["mode2_aligned"] = XY_aligned[X_pca.shape[0] :]

print("Write output to file")
adata.uns["method_id"] = "harmonic_alignment"
adata.write_h5ad(par["output"], compression = "gzip")

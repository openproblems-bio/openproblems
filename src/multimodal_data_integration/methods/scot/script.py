## VIASH START
par = {
  input = "output.h5ad",
  output = "output.scot.h5ad",
  n_svd = 100,
  balanced=False,
}
resources_dir = "../../utils/"
## VIASH END

print("Loading dependencies")
import scanpy as sc
import sklearn.decomposition
from SCOT import SCOT

# importing helper functions from common preprocessing.py file in resources dir
import sys
sys.path.append(resources_dir)
from preprocessing import log_cpm
from preprocessing import sqrt_cpm


print("Reading input h5ad file")
adata = sc.read_h5ad(par["input"])

print("Normalising mode 1")
sqrt_cpm(adata)

print("Normalising mode 2")
log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")


print("Performing PCA reduction")
n_svd = min([par["n_svd"], min(adata.X.shape) - 1, min(adata.obsm["mode2"].shape) - 1])
X_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.X)
Y_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.obsm["mode2"])

print("Initialize SCOT")
scot = SCOT(X_pca, Y_pca)

print("Call the unbalanced alignment")
# From https://github.com/rsinghlab/SCOT/blob/master/examples/unbalanced_GW_SNAREseq.ipynb # noqa: 501
X_new_unbal, y_new_unbal = scot.align(
    k=50, e=1e-3, rho=0.0005, normalize=True, balanced=par["balanced"]
)

print()
adata.obsm["aligned"] = X_new_unbal
adata.obsm["mode2_aligned"] = y_new_unbal

print("Write output to file")
adata.uns["method_id"] = "scot"
adata.write_h5ad(par["output"], compression = "gzip")

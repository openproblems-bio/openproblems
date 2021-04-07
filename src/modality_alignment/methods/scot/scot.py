## VIASH START
par = {
  input = "output.h5ad",
  output = "output.scot.h5ad",
  n_svd = 100,
  balanced=False,
}
resources_dir = "../../resources/utils.py"
## VIASH END

import scanpy as sc
import sklearn.decomposition
import sys
sys.path.append(resources_dir)

# importing helper functions from common preprocessing.py file in resources dir
from preprocessing import log_cpm
from preprocessing import sqrt_cpm

def _scot(adata, n_svd=100, balanced=False):
    from SCOT import SCOT

    # PCA reduction
    n_svd = min([n_svd, min(adata.X.shape) - 1, min(adata.obsm["mode2"].shape) - 1])
    X_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.X)
    Y_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.obsm["mode2"])

    # Initialize SCOT
    scot = SCOT(X_pca, Y_pca)

    # call the unbalanced alignment
    # From https://github.com/rsinghlab/SCOT/blob/master/examples/unbalanced_GW_SNAREseq.ipynb # noqa: 501
    X_new_unbal, y_new_unbal = scot.align(
        k=50, e=1e-3, rho=0.0005, normalize=True, balanced=balanced
    )
    adata.obsm["aligned"] = X_new_unbal
    adata.obsm["mode2_aligned"] = y_new_unbal

    return adata

if __name__ == "__main__":
    adata = sc.read_h5ad(par["input"])
    # Normalize mode1
    sqrt_cpm(adata)
    # Normalize mode2
    log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    # run scot
    _scot(adata, n_svd=par["n_svd"], balanced=par["balanced"])
    # Write output to file
    adata.write_h5ad(par["output"], compression=9)

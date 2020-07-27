from sklearn.decomposition import TruncatedSVD
import scipy.spatial
from ....tools.normalize import log_cpm


def procrustes(adata, n_svd=100):
    if min(adata.X.shape) <= n_svd:
        n_svd = min(adata.X.shape) - 1
    if min(adata.obsm["mode2"].shape) <= n_svd:
        n_svd = min(adata.obsm["mode2"].shape) - 1
    log_cpm(adata)
    X_pca = TruncatedSVD(n_svd).fit_transform(adata.X)
    Y_pca = TruncatedSVD(n_svd).fit_transform(adata.obsm["mode2"])
    X_proc, Y_proc, _ = scipy.spatial.procrustes(X_pca, Y_pca)
    adata.obsm["aligned"] = X_proc
    adata.obsm["mode2_aligned"] = Y_proc

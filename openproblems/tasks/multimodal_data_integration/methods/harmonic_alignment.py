import harmonicalignment
import numpy as np
from sklearn.decomposition import TruncatedSVD


def harmonic_alignment(adata, n_svd=100):
    if min(adata.X.shape) <= n_svd:
        n_svd = min(adata.X.shape) - 1
    if min(adata.obsm["mode2"].shape) <= n_svd:
        n_svd = min(adata.obsm["mode2"].shape) - 1
    X_pca = TruncatedSVD(n_svd).fit_transform(adata.X)
    Y_pca = TruncatedSVD(n_svd).fit_transform(adata.obsm["mode2"])
    ha_op = harmonicalignment.HarmonicAlignment(n_filters=8)
    ha_op.align(X_pca, Y_pca)
    XY_aligned = ha_op.diffusion_map()
    adata.obsm["aligned"] = XY_aligned[: X_pca.shape[0]]
    adata.obsm["mode2_aligned"] = XY_corrected[X_pca.shape[0] :]

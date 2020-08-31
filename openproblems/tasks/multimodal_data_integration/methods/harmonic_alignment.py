import numpy as np
from sklearn.decomposition import TruncatedSVD
from ....tools.normalize import sqrt_cpm, log_cpm, log_scran_pooling


def _harmonic_alignment(adata, n_svd=100):
    import harmonicalignment

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
    adata.obsm["mode2_aligned"] = XY_aligned[X_pca.shape[0] :]


def harmonic_alignment_sqrt_cpm(adata, n_svd=100):
    sqrt_cpm(adata)
    log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    _harmonic_alignment(adata, n_svd=n_svd)


def harmonic_alignment_log_scran_pooling(adata, n_svd=100):
    log_scran_pooling(adata)
    log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    _harmonic_alignment(adata, n_svd=n_svd)

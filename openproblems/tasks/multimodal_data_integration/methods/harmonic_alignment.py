import numpy as np
from sklearn.decomposition import TruncatedSVD
from ....tools.normalize import sqrt_cpm, log_cpm, log_scran_pooling
from ....tools.decorators import method
from ....tools.utils import check_version


def _harmonic_alignment(adata, n_svd=100, n_eigenvectors=100, n_pca_XY=100):
    import harmonicalignment

    n_svd = min([n_svd, min(adata.X.shape) - 1, min(adata.obsm["mode2"].shape) - 1])
    if adata.X.shape[0] <= n_eigenvectors:
        n_eigenvectors = None
    if adata.X.shape[0] <= n_pca_XY:
        n_pca_XY = None
    X_pca = TruncatedSVD(n_svd).fit_transform(adata.X)
    Y_pca = TruncatedSVD(n_svd).fit_transform(adata.obsm["mode2"])
    ha_op = harmonicalignment.HarmonicAlignment(
        n_filters=8, n_pca_XY=n_pca_XY, n_eigenvectors=n_eigenvectors
    )
    ha_op.align(X_pca, Y_pca)
    XY_aligned = ha_op.diffusion_map(n_eigenvectors=n_eigenvectors)
    adata.obsm["aligned"] = XY_aligned[: X_pca.shape[0]]
    adata.obsm["mode2_aligned"] = XY_aligned[X_pca.shape[0] :]


@method(
    method_name="Harmonic Alignment (sqrt CPM)",
    paper_name="Harmonic Alignment",
    paper_url="https://epubs.siam.org/doi/abs/10.1137/1.9781611976236.36",
    paper_year=2020,
    code_url="https://github.com/KrishnaswamyLab/harmonic-alignment",
    code_version=check_version("harmonicalignment"),
)
def harmonic_alignment_sqrt_cpm(adata, n_svd=100):
    sqrt_cpm(adata)
    log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    _harmonic_alignment(adata, n_svd=n_svd)


@method(
    method_name="Harmonic Alignment (log scran)",
    paper_name="Harmonic Alignment",
    paper_url="https://epubs.siam.org/doi/abs/10.1137/1.9781611976236.36",
    paper_year=2020,
    code_url="https://github.com/KrishnaswamyLab/harmonic-alignment",
    code_version=check_version("harmonicalignment"),
    image="openproblems-r-base",
)
def harmonic_alignment_log_scran_pooling(adata, n_svd=100):
    log_scran_pooling(adata)
    log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    _harmonic_alignment(adata, n_svd=n_svd)

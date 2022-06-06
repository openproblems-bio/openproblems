from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_scran_pooling
from ....tools.normalize import sqrt_cpm
from ....tools.utils import check_version

import sklearn.decomposition


def _harmonic_alignment(
    adata, test=False, n_svd=None, n_eigenvectors=None, n_pca_XY=None, n_filters=None
):
    import harmonicalignment

    if test:
        n_svd = n_svd or 20
        n_eigenvectors = n_eigenvectors or 20
        n_pca_XY = n_pca_XY or 20
        n_filters = n_filters or 4
    else:  # pragma: no cover
        n_svd = n_svd or 100
        n_eigenvectors = n_eigenvectors or 100
        n_pca_XY = n_pca_XY or 100
        n_filters = n_filters or 8
    n_svd = min([n_svd, min(adata.X.shape) - 1, min(adata.obsm["mode2"].shape) - 1])
    if adata.X.shape[0] <= n_eigenvectors:
        n_eigenvectors = None
    if adata.X.shape[0] <= n_pca_XY:
        n_pca_XY = None
    X_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.X)
    Y_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.obsm["mode2"])
    ha_op = harmonicalignment.HarmonicAlignment(
        n_filters=n_filters, n_pca_XY=n_pca_XY, n_eigenvectors=n_eigenvectors
    )
    ha_op.align(X_pca, Y_pca)
    XY_aligned = ha_op.diffusion_map(n_eigenvectors=n_eigenvectors)
    adata.obsm["aligned"] = XY_aligned[: X_pca.shape[0]]
    adata.obsm["mode2_aligned"] = XY_aligned[X_pca.shape[0] :]

    adata.uns["method_code_version"] = check_version("harmonicalignment")
    return adata


@method(
    method_name="Harmonic Alignment (sqrt CPM)",
    paper_name="Harmonic Alignment",
    paper_url="https://epubs.siam.org/doi/abs/10.1137/1.9781611976236.36",
    paper_year=2020,
    code_url="https://github.com/KrishnaswamyLab/harmonic-alignment",
    image="openproblems-python-extras",
)
def harmonic_alignment_sqrt_cpm(
    adata, test=False, n_svd=None, n_eigenvectors=None, n_pca_XY=None, n_filters=None
):
    adata = sqrt_cpm(adata)
    adata = log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    _harmonic_alignment(
        adata,
        test=test,
        n_svd=n_svd,
        n_eigenvectors=n_eigenvectors,
        n_pca_XY=n_pca_XY,
        n_filters=n_filters,
    )
    return adata


@method(
    method_name="Harmonic Alignment (log scran)",
    paper_name="Harmonic Alignment",
    paper_url="https://epubs.siam.org/doi/abs/10.1137/1.9781611976236.36",
    paper_year=2020,
    code_url="https://github.com/KrishnaswamyLab/harmonic-alignment",
    code_version=check_version("harmonicalignment"),
    image="openproblems-r-extras",
)
def harmonic_alignment_log_scran_pooling(
    adata, test=False, n_svd=None, n_eigenvectors=None, n_pca_XY=None, n_filters=None
):
    adata = log_scran_pooling(adata)
    adata = log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    _harmonic_alignment(
        adata,
        test=test,
        n_svd=n_svd,
        n_eigenvectors=n_eigenvectors,
        n_pca_XY=n_pca_XY,
        n_filters=n_filters,
    )
    return adata

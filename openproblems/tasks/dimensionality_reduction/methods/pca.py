from ....tools.decorators import method
from ....tools.normalize import log_cp10k
from ....tools.normalize import log_cp10k_hvg
from ....tools.utils import check_version

import functools

_pca_method = functools.partial(
    method,
    method_summary=(
        'PCA or "Principal Component Analysis" is a linear method that finds orthogonal'
        " directions in the data that capture the most variance. The first two"
        " principal components are chosen as the two-dimensional embedding. We select"
        " only the first two principal components as the two-dimensional embedding. PCA"
        " is calculated on the logCPM expression matrix with and without selecting 1000"
        " HVGs."
    ),
    paper_name="On lines and planes of closest fit to systems of points in space",
    paper_reference="pearson1901pca",
    paper_year=1901,
    code_url=(
        "https://scikit-learn.org/stable/modules/generated/"
        "sklearn.decomposition.PCA.html"
    ),
)


def _pca(adata, genes=None):
    import scanpy as sc

    if genes is not None:
        X = adata[:, genes].copy().X
    else:
        X = adata.X

    adata.obsm["X_emb"] = sc.tl.pca(X, n_comps=2, svd_solver="arpack")
    adata.uns["method_code_version"] = check_version("scikit-learn")
    return adata


@_pca_method(method_name="PCA (logCP10k)")
def pca_logCP10k(adata, test: bool = False):
    adata = log_cp10k(adata)
    return _pca(adata)


@_pca_method(method_name="PCA (logCP10k, 1kHVG)")
def pca_logCP10k_1kHVG(adata, test: bool = False):
    adata = log_cp10k_hvg(adata)
    return _pca(adata, genes=adata.var["highly_variable"])

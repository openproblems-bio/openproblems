from ....tools.decorators import method
from ....tools.normalize import log_cp10k
from ....tools.normalize import log_cp10k_hvg
from ....tools.utils import check_version

import functools

_pca_method = functools.partial(
    method,
    paper_name="On lines and planes of closest fit to systems of points in space",
    paper_reference="pearson1901pca",
    paper_year=1901,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.decomposition.PCA.html",
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


@_pca_method(method_name="Principle Component Analysis (PCA) (logCP10k)")
def pca_logCP10k(adata, test: bool = False):
    adata = log_cp10k(adata)
    return _pca(adata)


@_pca_method(method_name="Principle Component Analysis (PCA) (logCP10k, 1kHVG)")
def pca_logCP10k_1kHVG(adata, test: bool = False):
    adata = log_cp10k_hvg(adata)
    return _pca(adata, genes=adata.var["highly_variable"])

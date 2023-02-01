from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version

import functools

_pca_method = functools.partial(
    method,
    method_summary="TODO",
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


@_pca_method(method_name="Principle Component Analysis (PCA) (logCPM)")
def pca_logCPM(adata, test: bool = False):
    adata = log_cpm(adata)
    return _pca(adata)


@_pca_method(method_name="Principle Component Analysis (PCA) (logCPM, 1kHVG)")
def pca_logCPM_1kHVG(adata, test: bool = False):
    adata = log_cpm_hvg(adata)
    return _pca(adata, genes=adata.var["highly_variable"])

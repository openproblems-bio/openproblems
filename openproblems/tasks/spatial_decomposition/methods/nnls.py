from ....tools.decorators import method
from ....tools.utils import check_version
from .._utils import normalize_coefficients
from .._utils import obs_means
from .._utils import split_sc_and_sp

import numpy as np


@method(
    method_name="NNLS: Non-Negative Least Square",
    paper_name="Solving Least Squares Problems",
    paper_url="https://epubs.siam.org/doi/pdf/10.1137/1.9781611971217.bm",
    paper_year=1987,
    code_url="https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html",  # noqa: E501
)
def nnls_scipy(adata, test=False):
    from scipy.optimize import nnls
    from scipy.sparse import issparse

    adata_sc, adata = split_sc_and_sp(adata)
    adata_means = obs_means(adata_sc, "label")

    X = adata_means.X.T
    y = adata.X.T
    if issparse(X):
        X = X.toarray()
    if issparse(y):
        y = y.toarray()
    res = np.zeros((y.shape[1], X.shape[1]))  # (voxels,cells)
    for i in range(y.shape[1]):
        x, _ = nnls(X, y[:, i])
        res[i] = x

    res_prop = normalize_coefficients(res)

    adata.obsm["proportions_pred"] = res_prop
    adata.uns["method_code_version"] = check_version("scipy")

    return adata

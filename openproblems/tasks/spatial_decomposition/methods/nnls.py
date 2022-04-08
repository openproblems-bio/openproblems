from ....tools.decorators import method
from ....tools.utils import check_version
from .._utils import normalize_coefficients
from .._utils import obs_means

import numpy as np
import pandas as pd


@method(
    method_name="NNLS",
    paper_name="Solving Least Squares Problems",
    paper_url="https://epubs.siam.org/doi/pdf/10.1137/1.9781611971217.bm",
    paper_year=1987,
    code_url="https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html",  # noqa: E501
    code_version=check_version("scipy"),
)
def nnls_scipy(adata, test=False):
    from scipy.optimize import nnls
    from scipy.sparse import issparse

    adata_sc = adata.uns["sc_reference"].copy()
    labels = adata_sc.obs["label"].cat.categories
    adata_means = obs_means(adata_sc, "label")

    if issparse(adata.X):
        X = adata_means.X.T.toarray()
        y = adata.X.T.toarray()
    else:
        X = adata_means.X.T
        y = adata.X.T
    res = np.zeros((y.shape[1], X.shape[1]))  # (voxels,cells)
    for i in range(y.shape[1]):
        x, _ = nnls(X, y[:, i])
        res[i] = x

    res_prop = normalize_coefficients(res)

    adata.obsm["proportions_pred"] = pd.DataFrame(
        res_prop, columns=labels, index=adata.obs_names
    )

    return adata

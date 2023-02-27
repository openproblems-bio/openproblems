from ....tools.decorators import method
from ....tools.utils import check_version
from ..utils import normalize_coefficients
from ..utils import obs_means
from ..utils import split_sc_and_sp

import numpy as np


@method(
    method_name="Non-Negative Least Squares",
    method_summary=(
        "NNLS13 is a decomposition method based on Non-Negative Least Square Regression"
        " (NNLS). It was originally introduced by the method AutoGenes"
    ),
    paper_name=(
        "AutoGeneS: Automatic gene selection using multi-objective optimization for"
        " RNA-seq deconvolution"
    ),
    paper_reference="aliee2021autogenes",
    paper_year=2021,
    code_url=(
        "https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html"
    ),
)
def nnls_scipy(adata, test=False):
    from scipy.optimize import nnls
    from scipy.sparse import issparse

    adata_sc, adata = split_sc_and_sp(adata)
    adata_means = obs_means(adata_sc, "label")

    X = adata_means.X.T
    y = adata.X.T
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

from ....tools.decorators import method
from ....tools.utils import check_version
from ..utils import normalize_coefficients
from ..utils import obs_means
from ..utils import split_sc_and_sp
from typing import Optional

import numpy as np


@method(
    method_name="NuSVR",
    paper_name="Probabilistic Outputs for Support Vector Machines and Comparisons to Regularized Likelihood Methods",  # noqa: E501
    paper_url="http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.41.1639",
    paper_year=1999,
    code_url="https://scikit-learn.org/stable/modules/generated/sklearn.svm.NuSVR.html",
)
def nusvr_sklearn(adata, test: bool = False, max_iter: Optional[int] = None):
    from scipy.sparse import issparse
    from sklearn.svm import NuSVR

    if test:
        max_iter = max_iter or 100
    else:  # pragma: nocover
        max_iter = max_iter or 10000

    adata_sc, adata = split_sc_and_sp(adata)
    adata_means = obs_means(adata_sc, "label")

    X = adata_means.X.T
    y = adata.X.T
    if issparse(y):
        y = y.toarray()
    res = np.zeros((y.shape[1], X.shape[1]))  # (voxels,cells)
    for i in range(y.shape[1]):
        model = NuSVR(kernel="linear", max_iter=max_iter)
        model.fit(X, y[:, i])
        res[i] = model.coef_

    res_prop = normalize_coefficients(res)

    adata.obsm["proportions_pred"] = res_prop
    adata.uns["method_code_version"] = check_version("scikit-learn")
    return adata

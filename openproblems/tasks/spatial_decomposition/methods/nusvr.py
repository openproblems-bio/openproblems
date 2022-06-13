from ....tools.decorators import method
from ....tools.utils import check_version
from .._utils import normalize_coefficients
from .._utils import obs_means
from .._utils import split_sc_and_sp

import numpy as np
import pandas as pd


@method(
    method_name="NuSVR",
    paper_name="Probabilistic Outputs for Support Vector Machines and Comparisons to Regularized Likelihood Methods",  # noqa: E501
    paper_url="http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.41.1639",
    paper_year=1999,
    code_url="https://scikit-learn.org/stable/modules/generated/sklearn.svm.NuSVR.html",
    code_version=check_version("scikit-learn"),
)
def nusvr_sklearn(adata, test=False):
    from scipy.sparse import issparse
    from sklearn.svm import NuSVR

    adata_sc, adata = split_sc_and_sp(adata)
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
        model = NuSVR(kernel="linear")
        model.fit(X, y[:, i])
        res[i] = model.coef_

    res_prop = normalize_coefficients(res)

    adata.obsm["proportions_pred"] = pd.DataFrame(
        res_prop, columns=labels, index=adata.obs_names
    )
    return adata

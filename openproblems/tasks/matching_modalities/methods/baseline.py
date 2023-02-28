from ....tools.decorators import baseline_method
from ....tools.normalize import log_cp10k
from ....tools.utils import check_version

import numpy as np


@baseline_method(
    method_name="Random Features",
    method_summary=(
        "20-dimensional SVD is computed on the first modality, and is then randomly"
        " permuted twice, once for use as the output for each modality, producing"
        " random features with no correlation between modalities."
    ),
)
def random_features(adata, test=False, n_svd=20):
    import sklearn.decomposition

    n_svd = min([n_svd, min(adata.X.shape) - 1, min(adata.obsm["mode2"].shape) - 1])
    adata = log_cp10k(adata)
    X_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.X)
    adata.obsm["aligned"] = X_pca[np.random.permutation(np.arange(adata.shape[0]))]
    adata.obsm["mode2_aligned"] = X_pca[
        np.random.permutation(np.arange(adata.shape[0]))
    ]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="True Features",
    method_summary=(
        "20-dimensional SVD is computed on the first modality, and this same embedding"
        " is used as output for both modalities, producing perfectly aligned features"
        " from each modality."
    ),
)
def true_features(adata, test=False, n_svd=20):
    import sklearn.decomposition

    n_svd = min([n_svd, min(adata.X.shape) - 1, min(adata.obsm["mode2"].shape) - 1])
    adata = log_cp10k(adata)
    X_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.X)
    adata.obsm["aligned"] = X_pca
    adata.obsm["mode2_aligned"] = X_pca
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

from ....tools.decorators import metric
from scipy import sparse

import numpy as np
import scprep


def _square(X):
    if sparse.issparse(X):
        X.data = X.data**2
        return X
    else:
        return scprep.utils.toarray(X) ** 2


@metric(
    metric_name="Mean squared error",
    metric_summary=(
        "Mean squared error (MSE) is the average distance between each pair of matched "
        "observations of the same cell in the learned latent space. Lower is better."
    ),
    paper_reference="lance2022multimodal",
    maximize=False,
)
def mse(adata):
    X = scprep.utils.toarray(adata.obsm["aligned"])
    Y = scprep.utils.toarray(adata.obsm["mode2_aligned"])

    X_shuffled = X[np.random.permutation(np.arange(X.shape[0])), :]
    error_random = np.mean(np.sum(_square(X_shuffled - Y)))
    error_abs = np.mean(np.sum(_square(X - Y)))
    return error_abs / error_random

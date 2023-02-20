from ....tools.decorators import baseline_method
from ....tools.utils import check_version

import numpy as np


@baseline_method(
    method_name="Random Scores",
    method_summary="TODO",
)
def random_scores(adata, test=False):
    adata.obsm["gene_score"] = adata.X[
        np.random.permutation(np.arange(adata.X.shape[0]))
    ]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="True Scores",
    method_summary="TODO",
)
def true_scores(adata, test=False):
    adata.obsm["gene_score"] = adata.X
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

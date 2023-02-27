from ....tools.decorators import baseline_method
from ....tools.utils import check_version

import numpy as np


@baseline_method(
    method_name="Random Scores",
    method_summary=(
        "Random generation of gene scores by random permutation of gene expression"
        " values"
    ),
)
def random_scores(adata, test=False):
    adata.obsm["gene_score"] = adata.X[
        np.random.permutation(np.arange(adata.X.shape[0]))
    ]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="True Scores",
    method_summary="Perfect prediction of gene scores from gene expression values",
)
def true_scores(adata, test=False):
    adata.obsm["gene_score"] = adata.X
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

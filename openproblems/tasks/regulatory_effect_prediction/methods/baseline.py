from ....tools.decorators import method
from ....tools.utils import check_version

import numpy as np


@method(
    method_name="Random Scores",
    paper_name="Random Scores (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def random_scores(adata, test=False):
    adata.obsm["gene_score"] = adata.X[
        np.random.permutation(np.arange(adata.X.shape[0]))
    ]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@method(
    method_name="True Scores",
    paper_name="True Scores (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def true_scores(adata, test=False):
    adata.obsm["gene_score"] = adata.X
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

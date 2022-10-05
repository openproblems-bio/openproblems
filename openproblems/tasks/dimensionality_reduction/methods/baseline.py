from ....tools.decorators import method
from ....tools.utils import check_version

import numpy as np


@method(
    method_name="Random Features",
    paper_name="Random Features (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def random_features(adata, test=False):
    adata.obsm["X_emb"] = np.random.normal(0, 1, (adata.shape[0], 2))
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@method(
    method_name="All Features",
    paper_name="All Features (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def true_features(adata, test=False):
    adata.obsm["X_emb"] = adata.X
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

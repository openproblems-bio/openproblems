from .....tools.decorators import method
from .....tools.utils import check_version

import numpy as np


@method(
    method_name="No Integration",
    paper_name="No Integration (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def no_integration(adata, test=False):
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@method(
    method_name="Random Integration",
    paper_name="Random Integration (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def random_integration(adata, test=False):
    adata.X = adata.X[np.random.permutation(np.arange(adata.shape[0]))]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@method(
    method_name="Random Integration by Celltype",
    paper_name="Random Integration by Celltype (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def celltype_integration(adata, test=False):
    for batch_name in np.unique(adata.obs["batch"]):
        batch_idx = np.argwhere(adata.obs["batch"].to_numpy() == batch_name).flatten()
        adata.X[batch_idx] = adata.X[np.random.permutation(batch_idx)]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

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
    adata.obsm["aligned"] = adata.obsm["mode2"][
        np.random.permutation(np.arange(adata.shape[0]))
    ]
    adata.obsm["mode2_aligned"] = adata.obsm["mode2"]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@method(
    method_name="True Features",
    paper_name="True Features (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def true_features(adata, test=False):
    adata.obsm["aligned"] = adata.obsm["mode2"]
    adata.obsm["mode2_aligned"] = adata.obsm["mode2"]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

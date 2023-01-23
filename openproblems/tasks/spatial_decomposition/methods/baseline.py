from ....tools.decorators import method
from ....tools.utils import check_version
from ..utils import split_sc_and_sp

import numpy as np


@method(
    method_name="Random Proportions",
    paper_name="Random Proportions (baseline)",
    paper_reference="openproblems",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def random_proportions(adata, test=False):
    adata_sc, adata = split_sc_and_sp(adata)
    label_distribution = adata_sc.obs["label"].value_counts()
    adata.obsm["proportions_pred"] = np.random.dirichlet(
        label_distribution,
        size=adata.shape[0],
    )

    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@method(
    method_name="True Proportions",
    paper_name="True Proportions (baseline)",
    paper_reference="openproblems",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def true_proportions(adata, test=False):
    _, adata = split_sc_and_sp(adata)
    adata.obsm["proportions_pred"] = adata.obsm["proportions_true"]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

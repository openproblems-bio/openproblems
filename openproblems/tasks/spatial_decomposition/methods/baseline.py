from ....tools.decorators import baseline_method
from ....tools.utils import check_version
from ..utils import split_sc_and_sp

import numpy as np


@baseline_method(
    method_name="Random Proportions",
    method_summary=(
        "Random assignment of predicted celltype proportions from a Dirichlet"
        " distribution."
    ),
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


@baseline_method(
    method_name="True Proportions",
    method_summary=(
        "Perfect assignment of predicted celltype proportions from the ground truth."
    ),
)
def true_proportions(adata, test=False):
    _, adata = split_sc_and_sp(adata)
    adata.obsm["proportions_pred"] = adata.obsm["proportions_true"]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

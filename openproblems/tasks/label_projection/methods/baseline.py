from ....tools.decorators import method
from ....tools.utils import check_version

import numpy as np


@method(
    method_name="Majority Vote",
    paper_name="Majority Vote (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def majority_vote(adata, test=False):
    majority = adata.obs.labels[adata.obs.is_train].value_counts().index[0]
    adata.obs["labels_pred"] = np.nan
    adata.obs.loc[~adata.obs.is_train, "labels_pred"] = majority

    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@method(
    method_name="Random Labels",
    paper_name="Random Labels (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def random_labels(adata, test=False):
    label_distribution = adata.obs.labels[adata.obs.is_train].value_counts()
    label_distribution = label_distribution / label_distribution.sum()
    adata.obs["labels_pred"] = np.nan
    adata.obs.loc[~adata.obs.is_train, "labels_pred"] = np.random.choice(
        label_distribution.index,
        size=(~adata.obs.is_train).sum(),
        replace=True,
        p=label_distribution,
    )

    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@method(
    method_name="True Labels",
    paper_name="True Labels (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def true_labels(adata, test=False):
    adata.obs["labels_pred"] = adata.obs["labels"]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

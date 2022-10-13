from ....tools.decorators import method
from ....tools.utils import check_version
from ..utils import split_sc_and_sp

import numpy as np


@method(
    method_name="Random assignment (baseline)",
    paper_name="Open Problems for Single Cell Analysis",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
)
def random_proportion_assignment(adata, test=False):
    _, adata = split_sc_and_sp(adata)
    n_types = adata.obsm["proportions_true"].shape[1]
    props = np.random.dirichlet(
        np.ones(n_types),
        size=adata.shape[0],
    )

    adata.obsm["proportions_pred"] = props
    adata.uns["method_code_version"] = check_version("openproblems")

    return adata

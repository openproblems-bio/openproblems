from ....tools.decorators import method

import numpy as np
import pandas as pd


@method(
    method_name="Random assignment of porportion values. For baseline reference.",
    paper_name="NA",
    paper_url="NA",
    paper_year="NA",
    code_url="NA",
    code_version="NA",
)
def random_proportion_assignment(adata, test=False):
    n_types = adata.obsm["proportions_true"].shape[1]
    props = np.random.dirichlet(
        np.ones(n_types),
        size=adata.shape[0],
    )
    props = pd.DataFrame(
        props,
        columns=adata.obsm["proportions_true"].columns,
        index=adata.obs.index,
    )

    adata["proportions_pred"] = props

    return adata

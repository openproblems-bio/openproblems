from ....tools.decorators import method

import numpy as np


@method(
    method_name="Trivial test method",
    paper_name="N/A",
    paper_url="https://github.com/dburkhardt",
    paper_year=2021,
    code_url="https://github.com/dburkhardt",
    code_version="1.0.2",
)
def trivial_method(adata):
    adata.obs["labels_pred"] = np.full(shape=adata.n_obs, fill_value="unknown")
    return adata

from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.utils import check_r_version

import functools

_liana = r_function("liana.R", args="sce, method")

_sca_method = functools.partial(
    method,
    paper_name="Predicting cell-to-cell communication networks using NATMI",
    paper_url="https://www.nature.com/articles/s41467-020-18873-z",
    paper_year=2021,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


@_sca_method(
    method_name="NATMI",
)
def sca(adata):
    # log-normalize
    adata = log_cpm(adata)
    adata.layers['logcounts'] = adata.layers['log_cpm']
    del adata.layers['log_cpm']

    # Run LIANA
    liana_res = _liana(adata, method="natmi")
    # Format results
    liana_res["score"] = liana_res["edge_specificity"]
    liana_res.sort_values("score", ascending=False, inplace=True)
    adata.uns["ccc"] = liana_res

    adata.uns["method_code_version"] = check_r_version("liana")

    return adata

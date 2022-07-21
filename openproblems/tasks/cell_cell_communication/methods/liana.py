from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.utils import check_r_version

import functools

_r_liana = r_function("liana.R", args="sce, ...")

_liana_method = functools.partial(
    method,
    paper_name="Comparison of methods and resources for cell-cell "
    "communication inference from single-cell RNA-Seq data",
    paper_url="https://www.nature.com/articles/s41467-022-30755-0",
    paper_year=2022,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


@_liana_method(
    method_name="LIANA",
)
def liana(adata, score_col="aggregate_rank", asc=True, **kwargs):
    # log-normalize
    adata = log_cpm(adata)
    adata.layers["logcounts"] = adata.layers["log_cpm"]
    del adata.layers["log_cpm"]

    # Run LIANA
    liana_res = _r_liana(adata, **kwargs)

    # Format results
    liana_res["score"] = liana_res[score_col]
    liana_res.sort_values("score", ascending=asc, inplace=True)
    adata.uns["ccc"] = liana_res

    adata.uns["method_code_version"] = check_r_version("liana")

    return adata

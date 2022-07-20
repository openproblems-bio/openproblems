from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.utils import check_r_version

import functools

_liana = r_function("liana.R", args="sce, method")

_logfc_method = functools.partial(
    method,
    paper_name="Comparison of methods and resources for cell-cell "
               "communication inference from single-cell RNA-Seq data",
    paper_url="https://www.nature.com/articles/s41467-022-30755-0",
    paper_year=2022,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


@_logfc_method(
    method_name="Mean log2FC",
)
def logfc(adata):
    # log-normalize
    adata = log_cpm(adata)
    adata.layers['logcounts'] = adata.layers['log_cpm']
    del adata.layers['log_cpm']

    # Run LIANA
    liana_res = _liana(adata, method="logfc")
    # Format results
    liana_res["score"] = liana_res["logfc_comb"]
    liana_res.sort_values("score", ascending=False, inplace=True)
    adata.uns["ccc"] = liana_res

    adata.uns["method_code_version"] = check_r_version("liana")

    return adata

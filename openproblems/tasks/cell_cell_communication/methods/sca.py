from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.utils import check_r_version

import functools

_liana = r_function("liana.R", args="sce, method")

_sca_method = functools.partial(
    method,
    paper_name="SingleCellSignalR: inference of intercellular networks "
               "from single-cell transcriptomics.",
    paper_url="https://academic.oup.com/nar/article/48/10/e55/5810485",
    paper_year=2021,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


@_sca_method(
    method_name="SingleCellSignalR",
)
def sca(adata):
    # log-normalize
    adata = log_cpm(adata)
    adata.layers['logcounts'] = adata.layers['log_cpm']
    del adata.layers['log_cpm']

    # Run LIANA
    liana_res = _liana(adata, method="sca")
    # Format results
    liana_res["score"] = liana_res["LRscore"]
    liana_res.sort_values("score", ascending=False, inplace=True)
    adata.uns["ccc"] = liana_res

    adata.uns["method_code_version"] = check_r_version("liana")

    return adata

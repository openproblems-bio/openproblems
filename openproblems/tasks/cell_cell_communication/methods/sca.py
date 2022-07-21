from ....tools.decorators import method
from .liana import liana

import functools

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
    return liana(adata, method="sca", score_col="LRscore", asc=False)

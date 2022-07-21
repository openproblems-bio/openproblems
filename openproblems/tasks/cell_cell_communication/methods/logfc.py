from ....tools.decorators import method
from .liana import liana

import functools

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
def logfc(adata, test=False):
    return liana(adata, method="logfc", score_col="logfc_comb", asc=False)

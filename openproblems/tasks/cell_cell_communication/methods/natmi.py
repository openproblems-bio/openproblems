from ....tools.decorators import method
from .liana import liana

import functools

_natmi_method = functools.partial(
    method,
    paper_name="Predicting cell-to-cell communication networks using NATMI.",
    paper_url="https://www.nature.com/articles/s41467-020-18873-z",
    paper_year=2021,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


@_natmi_method(
    method_name="NATMI",
)
def natmi(adata):
    return liana(adata, method="natmi", score_col="edge_specificity", asc=False)

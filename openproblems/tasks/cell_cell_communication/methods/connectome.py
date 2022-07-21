from ....tools.decorators import method
from .liana import liana

import functools

_connectome_method = functools.partial(
    method,
    paper_name="Computation and visualization of cell–cell signaling topologies in single-cell systems data using Connectome",
    paper_url="https://www.nature.com/articles/s41598-022-07959-x",
    paper_year=2022,
    code_url="https://github.com/saezlab/liana",
    image="openproblems-r-extras",
)


@_connectome_method(
    method_name="Connectome",
)
def connectome(adata):
    return liana(adata, method="connectome", score_col="weight_sc", asc=False)

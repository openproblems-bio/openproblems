from .....tools.decorators import method
from .....tools.utils import check_version
from typing import Optional

import functools

_scalex_method = functools.partial(
    method,
    paper_name="Online single-cell data integration through projecting heterogeneous "
    "datasets into a common cell-embedding space",
    paper_url="https://doi.org/10.1038/s41467-022-33758-z",
    paper_year=2022,
    code_url="https://github.com/jsxlei/SCALEX",
    image="openproblems-python-extras",
)


def _scalex(
    adata,
    test: bool = False,
    max_iteration: Optional[int] = None,
    min_features: Optional[int] = None,
    min_cells: Optional[int] = None,
    compute_neighbors: bool = False,
    compute_features: bool = False,
):
    import scalex
    import scanpy as sc

    if test:
        max_iteration = max_iteration or 2
    else:  # pragma: nocover
        max_iteration = max_iteration or 30000

    if test or compute_features:
        min_features = min_features or 1
    else:  # pragma: nocover
        min_features = min_features or 600

    min_cells = min_cells or 1

    adata = scalex.SCALEX(
        adata,
        batch_key="batch",
        ignore_umap=True,
        impute=compute_features,
        max_iteration=max_iteration,
        min_features=min_features,
        min_cells=min_cells,
    )
    adata.obsm["X_emb"] = adata.obsm["latent"]
    if compute_features:
        adata.X = adata.layers["impute"]
    if compute_neighbors:
        sc.pp.neighbors(adata, use_rep="X_emb")
    adata.uns["method_code_version"] = check_version("scalex")
    return adata


@_scalex_method(method_name="SCALEX")
def scalex(adata, test: bool = False, max_iteration: Optional[int] = None):
    return _scalex(
        adata,
        test=test,
        max_iteration=max_iteration,
        compute_neighbors=True,
        compute_features=False,
    )

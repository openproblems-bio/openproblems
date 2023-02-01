from .....tools.decorators import method
from .....tools.utils import check_version
from typing import Optional

import functools

_scalex_method = functools.partial(
    method,
    method_summary="TODO",
    paper_name=(
        "Online single-cell data integration through projecting heterogeneous "
        "datasets into a common cell-embedding space"
    ),
    paper_reference="xiong2021online",
    paper_year=2022,
    code_url="https://github.com/jsxlei/SCALEX",
    image="openproblems-python-pytorch",
)


def _scalex(
    adata,
    test: bool = False,
    n_top_features: int = 0,
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
        impute=adata.obs["batch"].cat.categories[0] if compute_features else False,
        max_iteration=max_iteration,
        min_features=min_features,
        min_cells=min_cells,
        n_top_features=n_top_features,
        outdir=None,
    )
    adata.obsm["X_emb"] = adata.obsm["latent"]
    if compute_features:
        adata.X = adata.layers["impute"]
    if compute_neighbors:
        sc.pp.neighbors(adata, use_rep="X_emb")
    adata.uns["method_code_version"] = check_version("scalex")
    return adata


@_scalex_method(method_name="SCALEX (full)")
def scalex_full(adata, test: bool = False, max_iteration: Optional[int] = None):
    return _scalex(
        adata,
        test=test,
        max_iteration=max_iteration,
        compute_neighbors=True,
        compute_features=False,
        n_top_features=0,
    )


@_scalex_method(method_name="SCALEX (hvg)")
def scalex_hvg(adata, test: bool = False, max_iteration: Optional[int] = None):
    return _scalex(
        adata,
        test=test,
        max_iteration=max_iteration,
        compute_neighbors=True,
        compute_features=False,
        n_top_features=2000,
    )

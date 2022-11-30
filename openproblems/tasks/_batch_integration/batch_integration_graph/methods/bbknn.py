from .....tools.decorators import method
from .....tools.utils import check_version
from typing import Optional

import functools

_bbknn_method = functools.partial(
    method,
    paper_name="BBKNN: fast batch alignment of single cell transcriptomes",
    paper_reference="polanski2020bbknn",
    paper_year=2020,
    code_url="https://github.com/Teichlab/bbknn",
    image="openproblems-python-batch-integration",  # only if required
)


def _run_bbknn(
    adata,
    batch: str,
    test: bool = False,
    n_pca: Optional[int] = None,
    annoy_n_trees: Optional[int] = None,
    neighbors_within_batch: Optional[int] = None,
):
    from scanpy.preprocessing import pca

    import bbknn

    kwargs = dict(batch_key=batch, copy=True)

    if test:
        n_pca = n_pca or 10
        kwargs["annoy_n_trees"] = annoy_n_trees or 1
        kwargs["neighbors_within_batch"] = neighbors_within_batch or 3
    else:  # pragma: no cover
        n_pca = n_pca or 50
        kwargs["annoy_n_trees"] = annoy_n_trees or 10
        default_neighbors_within_batch = 25 if adata.n_obs >= 1e5 else 3
        kwargs["neighbors_within_batch"] = (
            neighbors_within_batch or default_neighbors_within_batch
        )

    pca(adata, n_comps=n_pca, svd_solver="arpack")
    adata = bbknn.bbknn(adata, **kwargs)

    adata.uns["method_code_version"] = check_version("bbknn")
    return adata


@_bbknn_method(method_name="BBKNN (full/unscaled)")
def bbknn_full_unscaled(
    adata,
    test: bool = False,
    n_pca: Optional[int] = None,
    annoy_n_trees: Optional[int] = None,
    neighbors_within_batch: Optional[int] = None,
):
    adata = _run_bbknn(
        adata,
        "batch",
        test=test,
        n_pca=n_pca,
        annoy_n_trees=annoy_n_trees,
        neighbors_within_batch=neighbors_within_batch,
    )
    # Complete the result in-place
    return adata


@_bbknn_method(method_name="BBKNN (hvg/unscaled)")
def bbknn_hvg_unscaled(
    adata,
    test: bool = False,
    n_pca: Optional[int] = None,
    annoy_n_trees: Optional[int] = None,
    neighbors_within_batch: Optional[int] = None,
):
    from ._utils import hvg_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = _run_bbknn(
        adata,
        "batch",
        test=test,
        n_pca=n_pca,
        annoy_n_trees=annoy_n_trees,
        neighbors_within_batch=neighbors_within_batch,
    )
    # Complete the result in-place
    return adata


@_bbknn_method(method_name="BBKNN (hvg/scaled)")
def bbknn_hvg_scaled(
    adata,
    test: bool = False,
    n_pca: Optional[int] = None,
    annoy_n_trees: Optional[int] = None,
    neighbors_within_batch: Optional[int] = None,
):
    from ._utils import hvg_batch
    from ._utils import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = _run_bbknn(
        adata,
        "batch",
        test=test,
        n_pca=n_pca,
        annoy_n_trees=annoy_n_trees,
        neighbors_within_batch=neighbors_within_batch,
    )
    # Complete the result in-place
    return adata


@_bbknn_method(method_name="BBKNN (full/scaled)")
def bbknn_full_scaled(
    adata,
    test: bool = False,
    n_pca: Optional[int] = None,
    annoy_n_trees: Optional[int] = None,
    neighbors_within_batch: Optional[int] = None,
):
    from ._utils import scale_batch

    adata = scale_batch(adata, "batch")
    adata = _run_bbknn(
        adata,
        "batch",
        test=test,
        n_pca=n_pca,
        annoy_n_trees=annoy_n_trees,
        neighbors_within_batch=neighbors_within_batch,
    )
    # Complete the result in-place
    return adata

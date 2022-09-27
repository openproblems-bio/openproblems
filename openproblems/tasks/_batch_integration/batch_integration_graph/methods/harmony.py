from .....tools.decorators import method
from .....tools.utils import check_version
from typing import Optional

import functools

_harmony_method = functools.partial(
    method,
    paper_name="Fast, sensitive and accurate integration "
    "of single-cell data with Harmony",
    paper_url="https://www.nature.com/articles/s41592-019-0619-0",
    paper_year=2019,
    code_url="https://github.com/lilab-bcb/harmony-pytorch",
    image="openproblems-r-pytorch",
)


def _harmony(
    adata,
    batch: str,
    test: bool = False,
    n_pca: Optional[int] = None,
    max_iter_harmony: Optional[int] = None,
    max_iter_cluster: Optional[int] = None,
):
    from harmony import harmonize
    from openproblems.tools.normalize import log_scran_pooling

    import scanpy as sc

    if test:
        n_pca = n_pca or 10
        max_iter_harmony = max_iter_harmony or 3
        max_iter_cluster = max_iter_cluster or 3
    else:  # pragma: nocover
        n_pca = n_pca or 50
        max_iter_harmony = max_iter_harmony or 10
        max_iter_cluster = max_iter_cluster or 20

    sc.pp.pca(adata, n_comps=n_pca)
    adata.obsm["X_emb"] = harmonize(
        adata.obsm["X_pca"],
        adata.obs,
        batch_key=batch,
        max_iter_harmony=max_iter_harmony,
        max_iter_clustering=max_iter_cluster,
    )
    sc.pp.neighbors(adata, use_rep="X_emb")

    adata.uns["method_code_version"] = check_version("harmony-pytorch")
    return adata


@_harmony_method(method_name="Harmony (full/unscaled)")
def harmony_full_unscaled(
    adata,
    test: bool = False,
    n_pca: Optional[int] = None,
    max_iter_harmony: Optional[int] = None,
    max_iter_cluster: Optional[int] = None,
):
    return _harmony(
        adata,
        "batch",
        test=test,
        n_pca=n_pca,
        max_iter_harmony=max_iter_harmony,
        max_iter_cluster=max_iter_cluster,
    )


@_harmony_method(method_name="Harmony (hvg/unscaled)")
def harmony_hvg_unscaled(
    adata,
    test: bool = False,
    n_pca: Optional[int] = None,
    max_iter_harmony: Optional[int] = None,
    max_iter_cluster: Optional[int] = None,
):
    from ._utils import hvg_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    return _harmony(
        adata,
        "batch",
        test=test,
        n_pca=n_pca,
        max_iter_harmony=max_iter_harmony,
        max_iter_cluster=max_iter_cluster,
    )


@_harmony_method(method_name="Harmony (hvg/scaled)")
def harmony_hvg_scaled(
    adata,
    test: bool = False,
    n_pca: Optional[int] = None,
    max_iter_harmony: Optional[int] = None,
    max_iter_cluster: Optional[int] = None,
):
    from ._utils import hvg_batch
    from ._utils import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    return _harmony(
        adata,
        "batch",
        test=test,
        n_pca=n_pca,
        max_iter_harmony=max_iter_harmony,
        max_iter_cluster=max_iter_cluster,
    )


@_harmony_method(method_name="Harmony (full/scaled)")
def harmony_full_scaled(
    adata,
    test: bool = False,
    n_pca: Optional[int] = None,
    max_iter_harmony: Optional[int] = None,
    max_iter_cluster: Optional[int] = None,
):
    from ._utils import scale_batch

    adata = scale_batch(adata, "batch")
    return _harmony(
        adata,
        "batch",
        test=test,
        n_pca=n_pca,
        max_iter_harmony=max_iter_harmony,
        max_iter_cluster=max_iter_cluster,
    )

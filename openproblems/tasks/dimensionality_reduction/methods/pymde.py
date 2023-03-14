from ....tools.decorators import method
from ....tools.normalize import log_cp10k
from ....tools.normalize import log_cp10k_hvg
from ....tools.utils import check_version
from typing import Optional

import functools
import scanpy as sc

_pymde_method = functools.partial(
    method,
    method_summary=(
        "PyMDE is a Python implementation of minimum-distortion embedding. It is a"
        " non-linear method that preserves distances between cells or neighborhoods in"
        " the high-dimensional space. It is computed with options to preserve distances"
        " between cells or neighbourhoods and with the logCPM matrix with and without"
        " HVG selection as input."
    ),
    paper_name="Minimum-Distortion Embedding",
    paper_reference="agrawal2021mde",
    paper_year=2021,
    code_url="https://pymde.org/",
    image="openproblems-python-pytorch",
)


def _pymde(
    adata,
    method: str = "neighbors",
    genes=None,
    n_pca: Optional[int] = None,
    test: bool = False,
    max_iter: Optional[int] = None,
    memory_size: Optional[int] = None,
):
    import pymde

    if genes is not None:
        adata_input = adata[:, genes].copy()
    else:
        adata_input = adata

    embed_kwargs = {}
    if test:
        n_pca = n_pca or 20
        embed_kwargs["max_iter"] = max_iter or 20
        embed_kwargs["memory_size"] = memory_size or 2
    else:  # pragma: nocover
        n_pca = n_pca or 100
        if max_iter is not None:
            embed_kwargs["max_iter"] = max_iter
        if memory_size is not None:
            embed_kwargs["memory_size"] = memory_size
    sc.tl.pca(adata_input, n_comps=n_pca, svd_solver="arpack")
    X = adata_input.obsm["X_pca"]
    if method == "neighbors":
        mde_fn = pymde.preserve_neighbors
    elif method == "distances":
        mde_fn = pymde.preserve_distances
    else:
        raise NotImplementedError
    adata.obsm["X_emb"] = (
        mde_fn(X, embedding_dim=2, verbose=True)
        .embed(**embed_kwargs, verbose=True)
        .detach()
        .numpy()
    )
    adata.uns["method_code_version"] = check_version("pymde")
    return adata


@_pymde_method(
    method_name="PyMDE Preserve Neighbors (logCP10k)",
)
def pymde_neighbors_log_cp10k(
    adata,
    test: bool = False,
    max_iter: Optional[int] = None,
    memory_size: Optional[int] = None,
):
    adata = log_cp10k(adata)
    return _pymde(
        adata, method="neighbors", test=test, max_iter=max_iter, memory_size=memory_size
    )


@_pymde_method(
    method_name="PyMDE Preserve Neighbors (logCP10k, 1kHVG)",
)
def pymde_neighbors_log_cp10k_hvg(
    adata,
    test: bool = False,
    max_iter: Optional[int] = None,
    memory_size: Optional[int] = None,
):
    adata = log_cp10k_hvg(adata)
    return _pymde(
        adata,
        method="neighbors",
        genes=adata.var["highly_variable"],
        test=test,
        max_iter=max_iter,
        memory_size=memory_size,
    )


@_pymde_method(
    method_name="PyMDE Preserve Distances (logCP10k)",
)
def pymde_distances_log_cp10k(
    adata,
    test: bool = False,
    max_iter: Optional[int] = None,
    memory_size: Optional[int] = None,
):
    adata = log_cp10k(adata)
    return _pymde(
        adata, method="distances", test=test, max_iter=max_iter, memory_size=memory_size
    )


@_pymde_method(
    method_name="PyMDE Preserve Distances (logCP10k, 1kHVG)",
)
def pymde_distances_log_cp10k_hvg(
    adata,
    test: bool = False,
    max_iter: Optional[int] = None,
    memory_size: Optional[int] = None,
):
    adata = log_cp10k_hvg(adata)
    return _pymde(
        adata,
        method="distances",
        genes=adata.var["highly_variable"],
        test=test,
        max_iter=max_iter,
        memory_size=memory_size,
    )

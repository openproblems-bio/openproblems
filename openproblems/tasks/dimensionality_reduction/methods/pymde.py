from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version
from typing import Optional

import scanpy as sc


def _pymde(
    adata,
    method: str = "neighbors",
    test: bool = False,
    max_iter: Optional[int] = None,
    memory_size: Optional[int] = None,
):
    import pymde

    embed_kwargs = {}
    if test:
        sc.tl.pca(adata, n_comps=20, svd_solver="arpack")
        X = adata.obsm["X_pca"]
        embed_kwargs["max_iter"] = max_iter or 20
        embed_kwargs["memory_size"] = memory_size or 2
    else:
        X = adata.X
        if max_iter is not None:
            embed_kwargs["max_iter"] = max_iter
        if memory_size is not None:
            embed_kwargs["memory_size"] = memory_size
    if method == "neighbors":
        mde_fn = pymde.preserve_neighbors
    elif method == "distances":
        mde_fn = pymde.preserve_distances
    else:
        raise NotImplementedError
    adata.obsm["X_emb"] = mde_fn(X, embedding_dim=2, verbose=True).embed(
        **embed_kwargs, verbose=True
    )
    adata.uns["method_code_version"] = check_version("pymde")
    return adata


@method(
    method_name="PyMDE Preserve Neighbors (logCPM)",
    paper_name="Minimum-Distortion Embedding",
    paper_url="http://dx.doi.org/10.1561/2200000090",
    paper_year=2021,
    code_url="https://pymde.org/",
    image="openproblems-python-pytorch",
)
def pymde_neighbors_log_cpm(
    adata,
    test: bool = False,
    max_iter: Optional[int] = None,
    memory_size: Optional[int] = None,
):
    adata = log_cpm(adata)
    return _pymde(
        adata, method="neighbors", test=test, max_iter=max_iter, memory_size=memory_size
    )


@method(
    method_name="PyMDE Preserve Neighbors (logCPM, 1kHVG)",
    paper_name="Minimum-Distortion Embedding",
    paper_url="http://dx.doi.org/10.1561/2200000090",
    paper_year=2021,
    code_url="https://pymde.org/",
    image="openproblems-python-pytorch",
)
def pymde_neighbors_log_cpm_hvg(
    adata,
    test: bool = False,
    max_iter: Optional[int] = None,
    memory_size: Optional[int] = None,
):
    adata = log_cpm_hvg(adata)
    return _pymde(
        adata, method="neighbors", test=test, max_iter=max_iter, memory_size=memory_size
    )


@method(
    method_name="PyMDE Preserve Distances (logCPM)",
    paper_name="Minimum-Distortion Embedding",
    paper_url="http://dx.doi.org/10.1561/2200000090",
    paper_year=2021,
    code_url="https://pymde.org/",
    image="openproblems-python-pytorch",
)
def pymde_distances_log_cpm(
    adata,
    test: bool = False,
    max_iter: Optional[int] = None,
    memory_size: Optional[int] = None,
):
    adata = log_cpm(adata)
    return _pymde(
        adata, method="distances", test=test, max_iter=max_iter, memory_size=memory_size
    )


@method(
    method_name="PyMDE Preserve Distances (logCPM, 1kHVG)",
    paper_name="Minimum-Distortion Embedding",
    paper_url="http://dx.doi.org/10.1561/2200000090",
    paper_year=2021,
    code_url="https://pymde.org/",
    image="openproblems-python-pytorch",
)
def pymde_distances_log_cpm_hvg(
    adata,
    test: bool = False,
    max_iter: Optional[int] = None,
    memory_size: Optional[int] = None,
):
    adata = log_cpm_hvg(adata)
    return _pymde(
        adata, method="distances", test=test, max_iter=max_iter, memory_size=memory_size
    )

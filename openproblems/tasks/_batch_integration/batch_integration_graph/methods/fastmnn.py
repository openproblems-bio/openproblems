from .....tools.conversion import r_function
from .....tools.decorators import method
from .....tools.utils import check_r_version
from typing import Optional

import functools

_fastmnn_method = functools.partial(
    method,
    paper_name="A description of the theory behind the fastMNN algorithm",
    paper_url="https://marionilab.github.io/FurtherMNN2018/theory/description.html",
    paper_year=2019,
    code_url="https://doi.org/doi:10.18129/B9.bioc.batchelor",
    image="openproblems-r-extras",
)

_r_fastmnn = r_function("fastmnn.R", args="sce, batch, k, n_pca, return_features=FALSE")


def _fastmnn(
    adata,
    batch,
    return_features: bool,
    test: bool = False,
    k: Optional[int] = None,
    n_pca: Optional[int] = None,
):
    from openproblems.tools.normalize import log_scran_pooling

    if test:
        k = k or 5
        n_pca = n_pca or 10
    else:  # pragma: nocover
        k = k or 20
        n_pca = n_pca or 50


    return _r_fastmnn(adata, batch, k, n_pca, return_features=return_features)


def _fastmnn_embed(
    adata,
    batch,
    test: bool = False,
    k: Optional[int] = None,
    n_pca: Optional[int] = None,
):
    from scanpy.preprocessing import neighbors

    adata.obsm["X_emb"] = _fastmnn(
        adata, batch, test=test, k=k, n_pca=n_pca, return_features=False
    )
    neighbors(adata, use_rep="X_emb")
    adata.uns["method_code_version"] = check_r_version("batchelor")
    return adata


def _fastmnn_feature(
    adata,
    batch,
    test: bool = False,
    k: Optional[int] = None,
    n_pca: Optional[int] = None,
):
    from scib.preprocessing import reduce_data

    adata.X = _fastmnn(
        adata, batch, test=test, k=k, n_pca=n_pca, return_features=True
    ).T
    reduce_data(adata, umap=False)
    adata.uns["method_code_version"] = check_r_version("batchelor")
    return adata


@_fastmnn_method(method_name="FastMNN embed (full/unscaled)")
def fastmnn_embed_full_unscaled(
    adata, test: bool = False, k: Optional[int] = None, n_pca: Optional[int] = None
):
    return _fastmnn_embed(adata, "batch", test=test, k=k, n_pca=n_pca)


@_fastmnn_method(method_name="FastMNN embed (hvg/unscaled)")
def fastmnn_embed_hvg_unscaled(
    adata, test: bool = False, k: Optional[int] = None, n_pca: Optional[int] = None
):
    from ._utils import hvg_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    return _fastmnn_embed(adata, "batch", test=test, k=k, n_pca=n_pca)


@_fastmnn_method(method_name="FastMNN embed (hvg/scaled)")
def fastmnn_embed_hvg_scaled(
    adata, test: bool = False, k: Optional[int] = None, n_pca: Optional[int] = None
):
    from ._utils import hvg_batch
    from ._utils import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    return _fastmnn_embed(adata, "batch", test=test, k=k, n_pca=n_pca)


@_fastmnn_method(method_name="FastMNN embed (full/scaled)")
def fastmnn_embed_full_scaled(
    adata, test: bool = False, k: Optional[int] = None, n_pca: Optional[int] = None
):
    from ._utils import scale_batch

    adata = scale_batch(adata, "batch")
    return _fastmnn_embed(adata, "batch", test=test, k=k, n_pca=n_pca)


@_fastmnn_method(method_name="FastMNN feature (full/unscaled)")
def fastmnn_feature_full_unscaled(
    adata, test: bool = False, k: Optional[int] = None, n_pca: Optional[int] = None
):
    return _fastmnn_feature(adata, "batch", test=test, k=k, n_pca=n_pca)


@_fastmnn_method(method_name="FastMNN feature (hvg/unscaled)")
def fastmnn_feature_hvg_unscaled(
    adata, test: bool = False, k: Optional[int] = None, n_pca: Optional[int] = None
):
    from ._utils import hvg_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    return _fastmnn_feature(adata, "batch", test=test, k=k, n_pca=n_pca)


@_fastmnn_method(method_name="FastMNN feature (hvg/scaled)")
def fastmnn_feature_hvg_scaled(
    adata, test: bool = False, k: Optional[int] = None, n_pca: Optional[int] = None
):
    from ._utils import hvg_batch
    from ._utils import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    return _fastmnn_feature(adata, "batch", test=test, k=k, n_pca=n_pca)


@_fastmnn_method(method_name="FastMNN feature (full/scaled)")
def fastmnn_feature_full_scaled(
    adata, test: bool = False, k: Optional[int] = None, n_pca: Optional[int] = None
):
    from ._utils import scale_batch

    adata = scale_batch(adata, "batch")
    return _fastmnn_feature(adata, "batch", test=test, k=k, n_pca=n_pca)

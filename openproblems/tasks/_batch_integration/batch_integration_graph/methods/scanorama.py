from .....tools.decorators import method
from .....tools.utils import check_version

import functools

_scanorama_method = functools.partial(
    method,
    paper_name="Efficient integration of heterogeneous single-cell "
    "transcriptomes using Scanorama",
    paper_url="https://www.nature.com/articles/s41587-019-0113-3",
    paper_year=2019,
    code_url="https://github.com/brianhie/scanorama",
    image="openproblems-python-batch-integration",
)


def _scanorama(adata, use_rep):
    from scib.integration import runScanorama
    from scib.preprocessing import reduce_data


    adata = runScanorama(adata, "batch")
    reduce_data(adata, umap=False, use_rep=use_rep)
    adata.uns["method_code_version"] = check_version("scanorama")
    return adata


def _scanorama_embed(adata):
    return _scanorama(adata, use_rep="X_emb")


def _scanorama_full(adata):
    return _scanorama(adata, use_rep="X_pca")


@_scanorama_method(method_name="Scanorama (full/unscaled)")
def scanorama_embed_full_unscaled(adata, test=False):
    return _scanorama_embed(adata)


@_scanorama_method(method_name="Scanorama (hvg/unscaled)")
def scanorama_embed_hvg_unscaled(adata, test=False):
    from ._utils import hvg_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    return _scanorama_embed(adata)


@_scanorama_method(method_name="Scanorama (hvg/scaled)")
def scanorama_embed_hvg_scaled(adata, test=False):
    from ._utils import hvg_batch
    from ._utils import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    return _scanorama_embed(adata)


@_scanorama_method(method_name="Scanorama (full/scaled)")
def scanorama_embed_full_scaled(adata, test=False):
    from ._utils import scale_batch

    adata = scale_batch(adata, "batch")
    return _scanorama_embed(adata)


@_scanorama_method(method_name="Scanorama gene output (full/unscaled)")
def scanorama_feature_full_unscaled(adata, test=False):
    return _scanorama_full(adata)


@_scanorama_method(method_name="Scanorama gene output (hvg/unscaled)")
def scanorama_feature_hvg_unscaled(adata, test=False):
    from ._utils import hvg_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    return _scanorama_full(adata)


@_scanorama_method(method_name="Scanorama gene output (hvg/scaled)")
def scanorama_feature_hvg_scaled(adata, test=False):
    from ._utils import hvg_batch
    from ._utils import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    return _scanorama_full(adata)


@_scanorama_method(method_name="Scanorama gene output (full/scaled)")
def scanorama_feature_full_scaled(adata, test=False):
    from ._utils import scale_batch

    adata = scale_batch(adata, "batch")
    return _scanorama_full(adata)

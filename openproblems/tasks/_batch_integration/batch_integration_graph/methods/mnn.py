from .....tools.decorators import method
from .....tools.utils import check_version

import functools

_mnn_method = functools.partial(
    method,
    paper_name="Batch effects in single-cell RNA-sequencing "
    "data are corrected by matching mutual nearest neighbors",
    paper_reference="haghverdi2018batch",
    paper_year=2018,
    code_url="https://github.com/chriscainx/mnnpy",
    image="openproblems-r-pytorch",
)


def _mnn(adata):
    from scib.integration import runMNN
    from scib.preprocessing import reduce_data

    # mnn clears adata.uns
    uns = adata.uns
    adata = runMNN(adata, "batch")
    adata.uns = uns
    reduce_data(adata, umap=False)
    adata.obsm["X_emb"] = adata.obsm["X_pca"]
    adata.uns["method_code_version"] = check_version("mnnpy")
    return adata


@_mnn_method(method_name="MNN (full/unscaled)")
def mnn_full_unscaled(adata, test=False):
    return _mnn(adata)


@_mnn_method(method_name="MNN (hvg/unscaled)")
def mnn_hvg_unscaled(adata, test=False):
    from ._utils import hvg_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    return _mnn(adata)


@_mnn_method(method_name="MNN (hvg/scaled)")
def mnn_hvg_scaled(adata, test=False):
    from ._utils import hvg_batch
    from ._utils import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    return _mnn(adata)


@_mnn_method(method_name="MNN (full/scaled)")
def mnn_full_scaled(adata, test=False):
    from ._utils import scale_batch

    adata = scale_batch(adata, "batch")
    return _mnn(adata)

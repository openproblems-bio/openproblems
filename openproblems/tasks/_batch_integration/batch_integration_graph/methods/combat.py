# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version


@method(
    method_name="MNN",
    paper_name="Sc"
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("mnnpy"),
    # image="openproblems-template-image" # only if required
)
def mnn_full_unscaled(adata):
    from scIB.integration import runMNN
    from scIB.preprocessing import hvg_batch
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch
    adata = runMNN(adata, "batch")
    reduce_data(adata)
    # Complete the result in-place
    return adata


def mnn_hvg_unscaled(adata):
    from scIB.integration import runMNN
    from scIB.preprocessing import hvg_batch
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = runMNN(adata, "batch")
    reduce_data(adata)
    return adata


def mnn_hvg_scaled(adata):
    from scIB.integration import runMNN
    from scIB.preprocessing import hvg_batch
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, 'batch')
    adata = runMNN(adata, "batch")
    reduce_data(adata)
    return adata


def mnn_full_scaled(adata):
    from scIB.integration import runMNN
    from scIB.preprocessing import hvg_batch
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch
    adata = scale_batch(adata, 'batch')
    adata = runMNN(adata, "batch")
    reduce_data(adata)
    return adata

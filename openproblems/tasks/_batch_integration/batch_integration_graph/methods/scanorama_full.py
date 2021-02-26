# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version
from scIB.integration import runScanorama
from scIB.preprocessing import hvg_batch
from scIB.preprocessing import reduce_data
from scIB.preprocessing import scale_batch


@method(
    method_name="Scanorama",
    paper_name="Sc"
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scanorama"),
    # image="openproblems-template-image" # only if required
)
def scanorama_feature_full_unscaled(adata):
    runScanorama(adata, "batch")
    reduce_data(adata)
    # Complete the result in-place
    return adata


def scanorama_feature_hvg_unscaled(adata):
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    runScanorama(adata, "batch")
    reduce_data(adata)
    return adata


def scanorama_feature_hvg_scaled(adata):
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, 'batch')
    runScanorama(adata, "batch")
    reduce_data(adata)
    return adata


def scanorama_feature_full_scaled(adata):
    adata = scale_batch(adata, 'batch')
    runScanorama(adata, "batch")
    reduce_data(adata)
    return adata

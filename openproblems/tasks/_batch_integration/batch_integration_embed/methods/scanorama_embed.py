# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version


@method(
    method_name="Scanorama",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scanorama"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanorama_embed_full_unscaled(adata):
    from scIB.integration import runScanorama

    adata = runScanorama(adata, "batch")
    # Complete the result in-place
    return adata


@method(
    method_name="Scanorama (hvg/unscaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scanorama"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanorama_embed_hvg_unscaled(adata):
    from _hvg import hvg_batch
    from scIB.integration import runScanorama

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = runScanorama(adata, "batch")
    return adata


@method(
    method_name="Scanorama (hvg/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scanorama"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanorama_embed_hvg_scaled(adata):
    from _hvg import hvg_batch
    from scIB.integration import runScanorama
    from scIB.preprocessing import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = runScanorama(adata, "batch")
    return adata


@method(
    method_name="Scanorama (full/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scanorama"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanorama_embed_full_scaled(adata):
    from scIB.integration import runScanorama
    from scIB.preprocessing import scale_batch

    adata = scale_batch(adata, "batch")
    adata = runScanorama(adata, "batch")
    return adata

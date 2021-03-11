# from ....tools.normalize import log_cpm
from .....tools.decorators import method

from .....tools.utils import check_version


@method(
    method_name="scanvi",
    paper_name="",
    paper_url="",
    paper_year=0,
    code_url="",
    code_version=check_version("scvi"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanvi_full_unscaled(adata):
    from scIB.integration import runScanvi
    from scIB.preprocessing import reduce_data

    adata.obs.rename(columns={"labels": "lab"})  # ugly fix for scvi conversion error
    adata = runScanvi(adata, "batch", "labels")
    reduce_data(adata, use_rep="X_emb")
    adata.obs.rename(columns={"lab": "labels"})  # ugly fix for scvi conversion error
    return adata


@method(
    method_name="Scanvi (hvg/unscaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scvi"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanvi_hvg_unscaled(adata):
    from scIB.integration import runScanvi
    from scIB.preprocessing import hvg_batch
    from scIB.preprocessing import reduce_data

    adata.obs.rename(columns={"labels": "lab"})  # ugly fix for scvi conversion error
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = runScanvi(adata, "batch", "lab")
    reduce_data(adata, use_rep="X_emb")
    adata.obs.rename(columns={"lab": "labels"})  # ugly fix for scvi conversion error
    return adata


@method(
    method_name="Scanvi (hvg/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scvi"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanvi_hvg_scaled(adata):
    from scIB.integration import runScanvi
    from scIB.preprocessing import hvg_batch
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch

    adata.obs.rename(columns={"labels": "lab"})  # ugly fix for scvi conversion error
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = runScanvi(adata, "batch", "lab")
    reduce_data(adata, use_rep="X_emb")
    adata.obs.rename(columns={"lab": "labels"})  # ugly fix for scvi conversion error
    return adata


@method(
    method_name="Scanvi (full/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scvi"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanvi_full_scaled(adata):
    from scIB.integration import runScanvi
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch

    adata.obs.rename(columns={"labels": "lab"})  # ugly fix for scvi conversion error
    adata = scale_batch(adata, "batch")
    adata = runScanvi(adata, "batch", "lab")
    reduce_data(adata, use_rep="X_emb")
    adata.obs.rename(columns={"lab": "labels"})  # ugly fix for scvi conversion error
    return adata

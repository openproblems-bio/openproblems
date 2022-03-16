# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version


@method(
    method_name="Saucie gene output",
    paper_name="Exploring single-cell data with deep multitasking neural networks",
    paper_url="https://www.nature.com/articles/s41592-019-0576-7",
    paper_year=2019,
    code_url="https://github.com/KrishnaswamyLab/SAUCIE",
    code_version=check_version("SAUCIE"),
    image="openproblems-python37-saucie",  # only if required
)
def saucie_feature_full_unscaled(adata):
    from scib.integration import runSaucie
    from scib.preprocessing import reduce_data

    adata = runSaucie(adata, "batch")
    adata.obs['batch']= adata.obs['batch'].astype('category')
    reduce_data(adata, umap=False)
    # Complete the result in-place
    return adata


@method(
    method_name="Saucie gene output (hvg/unscaled)",
    paper_name="Exploring single-cell data with deep multitasking neural networks",
    paper_url="https://www.nature.com/articles/s41592-019-0576-7",
    paper_year=2019,
    code_url="https://github.com/KrishnaswamyLab/SAUCIE",
    code_version=check_version("SAUCIE"),
    image="openproblems-python37-saucie",  # only if required
)
def saucie_feature_hvg_unscaled(adata):
    from ._utils import hvg_batch
    from scib.integration import runSaucie
    from scib.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata.obs['batch']= adata.obs['batch'].astype('category')
    adata = runSaucie(adata, "batch")
    reduce_data(adata, umap=False)
    return adata


@method(
    method_name="Saucie gene output (hvg/scaled)",
    paper_name="Exploring single-cell data with deep multitasking neural networks",
    paper_url="https://www.nature.com/articles/s41592-019-0576-7",
    paper_year=2019,
    code_url="https://github.com/KrishnaswamyLab/SAUCIE",
    code_version=check_version("SAUCIE"),
    image="openproblems-python37-saucie",  # only if required
)
def saucie_feature_hvg_scaled(adata):
    from ._utils import hvg_batch
    from scib.integration import runSaucie
    from scib.preprocessing import reduce_data
    from scib.preprocessing import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata.obs['batch']= adata.obs['batch'].astype('category')
    adata = runSaucie(adata, "batch")
    reduce_data(adata, umap=False)
    return adata


@method(
    method_name="Saucie gene output (full/scaled)",
    paper_name="Exploring single-cell data with deep multitasking neural networks",
    paper_url="https://www.nature.com/articles/s41592-019-0576-7",
    paper_year=2019,
    code_url="https://github.com/KrishnaswamyLab/SAUCIE",
    code_version=check_version("SAUCIE"),
    image="openproblems-python37-saucie",  # only if required
)
def saucie_feature_full_scaled(adata):
    from scib.integration import runSaucie
    from scib.preprocessing import reduce_data
    from scib.preprocessing import scale_batch

    adata = scale_batch(adata, "batch")
    adata.obs['batch']= adata.obs['batch'].astype('category')
    adata = runSaucie(adata, "batch")
    reduce_data(adata, umap=False)
    return adata

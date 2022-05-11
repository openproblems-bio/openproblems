# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version


@method(
    method_name="MNN",
    paper_name="""Batch effects in single-cell RNA-sequencing
               data are corrected by matching mutual nearest neighbors""",
    paper_url="https://www.nature.com/articles/nbt.4091",
    paper_year=2018,
    code_url="https://github.com/chriscainx/mnnpy",
    code_version=check_version("mnnpy"),
    image="openproblems-python-batch-integration",
)
def mnn_full_unscaled(adata, test=False):
    from scib.integration import runMNN
    from scib.preprocessing import reduce_data

    adata = runMNN(adata, "batch")
    reduce_data(adata, umap=False)
    # Complete the result in-place
    return adata


@method(
    method_name="MNN (hvg/unscaled)",
    paper_name="""Batch effects in single-cell RNA-sequencing
               data are corrected by matching mutual nearest neighbors""",
    paper_url="https://www.nature.com/articles/nbt.4091",
    paper_year=2018,
    code_url="https://github.com/chriscainx/mnnpy",
    code_version=check_version("mnnpy"),
    image="openproblems-python-batch-integration",
)
def mnn_hvg_unscaled(adata, test=False):
    from ._utils import hvg_batch
    from scib.integration import runMNN
    from scib.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = runMNN(adata, "batch")
    reduce_data(adata, umap=False)
    return adata


@method(
    method_name="MNN (hvg/scaled)",
    paper_name="""Batch effects in single-cell RNA-sequencing
               data are corrected by matching mutual nearest neighbors""",
    paper_url="https://www.nature.com/articles/nbt.4091",
    paper_year=2018,
    code_url="https://github.com/chriscainx/mnnpy",
    code_version=check_version("mnnpy"),
    image="openproblems-python-batch-integration",
)
def mnn_hvg_scaled(adata, test=False):
    from ._utils import hvg_batch
    from ._utils import scale_batch
    from scib.integration import runMNN
    from scib.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = runMNN(adata, "batch")
    reduce_data(adata, umap=False)
    return adata


@method(
    method_name="MNN (full/scaled)",
    paper_name="""Batch effects in single-cell RNA-sequencing
               data are corrected by matching mutual nearest neighbors""",
    paper_url="https://www.nature.com/articles/nbt.4091",
    paper_year=2018,
    code_url="https://github.com/chriscainx/mnnpy",
    code_version=check_version("mnnpy"),
    image="openproblems-python-batch-integration",
)
def mnn_full_scaled(adata, test=False):
    from ._utils import scale_batch
    from scib.integration import runMNN
    from scib.preprocessing import reduce_data

    adata = scale_batch(adata, "batch")
    adata = runMNN(adata, "batch")
    reduce_data(adata, umap=False)
    return adata

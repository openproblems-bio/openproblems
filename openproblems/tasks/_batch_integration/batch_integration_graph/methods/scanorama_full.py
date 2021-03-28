# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version


@method(
    method_name="Scanorama gene output",
    paper_name="""Efficient integration of heterogeneous single-cell
               transcriptomes using Scanorama""",
    paper_url="https://www.nature.com/articles/s41587-019-0113-3",
    paper_year=2019,
    code_url="https://github.com/brianhie/scanorama",
    code_version=check_version("scanorama"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanorama_feature_full_unscaled(adata):
    from scIB.integration import runScanorama
    from scIB.preprocessing import reduce_data

    adata = runScanorama(adata, "batch")
    reduce_data(adata, umap=False)
    # Complete the result in-place
    return adata


@method(
    method_name="Scanorama gene output (hvg/unscaled)",
    paper_name="""Efficient integration of heterogeneous single-cell
               transcriptomes using Scanorama""",
    paper_url="https://www.nature.com/articles/s41587-019-0113-3",
    paper_year=2019,
    code_url="https://github.com/brianhie/scanorama",
    code_version=check_version("scanorama"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanorama_feature_hvg_unscaled(adata):
    from ._hvg import hvg_batch
    from scIB.integration import runScanorama
    from scIB.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = runScanorama(adata, "batch")
    reduce_data(adata, umap=False)
    return adata


@method(
    method_name="Scanorama gene output (hvg/scaled)",
    paper_name="""Efficient integration of heterogeneous single-cell
               transcriptomes using Scanorama""",
    paper_url="https://www.nature.com/articles/s41587-019-0113-3",
    paper_year=2019,
    code_url="https://github.com/brianhie/scanorama",
    code_version=check_version("scanorama"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanorama_feature_hvg_scaled(adata):
    from ._hvg import hvg_batch
    from ._hvg import scale_batch
    from scIB.integration import runScanorama
    from scIB.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = runScanorama(adata, "batch")
    reduce_data(adata, umap=False)
    return adata


@method(
    method_name="Scanorama gene output (full/scaled)",
    paper_name="""Efficient integration of heterogeneous single-cell
               transcriptomes using Scanorama""",
    paper_url="https://www.nature.com/articles/s41587-019-0113-3",
    paper_year=2019,
    code_url="https://github.com/brianhie/scanorama",
    code_version=check_version("scanorama"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanorama_feature_full_scaled(adata):
    from ._hvg import scale_batch
    from scIB.integration import runScanorama
    from scIB.preprocessing import reduce_data

    adata = scale_batch(adata, "batch")
    adata = runScanorama(adata, "batch")
    reduce_data(adata, umap=False)
    return adata

# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version


@method(
    method_name="Scanorama",
    paper_name="""Efficient integration of heterogeneous single-cell
               transcriptomes using Scanorama""",
    paper_url="https://www.nature.com/articles/s41587-019-0113-3",
    paper_year=2019,
    code_url="https://github.com/brianhie/scanorama",
    code_version=check_version("scanorama"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanorama_embed_full_unscaled(adata, test=False):
    from scib.integration import runScanorama
    from scib.preprocessing import reduce_data

    adata = runScanorama(adata, "batch")
    reduce_data(adata, umap=False, use_rep="X_emb")
    # Complete the result in-place
    return adata


@method(
    method_name="Scanorama (hvg/unscaled)",
    paper_name="""Efficient integration of heterogeneous single-cell
               transcriptomes using Scanorama""",
    paper_url="https://www.nature.com/articles/s41587-019-0113-3",
    paper_year=2019,
    code_url="https://github.com/brianhie/scanorama",
    code_version=check_version("scanorama"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanorama_embed_hvg_unscaled(adata, test=False):
    from ._utils import hvg_batch
    from scib.integration import runScanorama
    from scib.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = runScanorama(adata, "batch")
    reduce_data(adata, umap=False, use_rep="X_emb")
    return adata


@method(
    method_name="Scanorama (hvg/scaled)",
    paper_name="""Efficient integration of heterogeneous single-cell
               transcriptomes using Scanorama""",
    paper_url="https://www.nature.com/articles/s41587-019-0113-3",
    paper_year=2019,
    code_url="https://github.com/brianhie/scanorama",
    code_version=check_version("scanorama"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanorama_embed_hvg_scaled(adata, test=False):
    from ._utils import hvg_batch
    from ._utils import scale_batch
    from scib.integration import runScanorama
    from scib.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = runScanorama(adata, "batch")
    reduce_data(adata, umap=False, use_rep="X_emb")
    return adata


@method(
    method_name="Scanorama (full/scaled)",
    paper_name="""Efficient integration of heterogeneous single-cell
               transcriptomes using Scanorama""",
    paper_url="https://www.nature.com/articles/s41587-019-0113-3",
    paper_year=2019,
    code_url="https://github.com/brianhie/scanorama",
    code_version=check_version("scanorama"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanorama_embed_full_scaled(adata, test=False):
    from ._utils import scale_batch
    from scib.integration import runScanorama
    from scib.preprocessing import reduce_data

    adata = scale_batch(adata, "batch")
    adata = runScanorama(adata, "batch")
    reduce_data(adata, umap=False, use_rep="X_emb")
    return adata

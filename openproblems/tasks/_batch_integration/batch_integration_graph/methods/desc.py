# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version


@method(
    method_name="DESC",
    paper_name="""Deep learning enables accurate clustering with
     batch effect removal in single-cell RNA-seq analysis""",
    paper_url="https://www.nature.com/articles/s41467-020-15851-3",
    paper_year=2020,
    code_url="https://github.com/eleozzr/desc",
    code_version=check_version("desc"),
    image="openproblems-python37-scgen",  # only if required
)
def desc_full_unscaled(adata):
    from scanpy.preprocessing import neighbors
    from scib.integration import runDESC

    adata.obsm["X_emb"] = runDESC(adata, "batch").obsm["X_emb"]
    neighbors(adata, use_rep="X_emb")
    # Complete the result in-place
    return adata


@method(
    method_name="DESC (hvg/unscaled)",
    paper_name="""Deep learning enables accurate clustering with
     batch effect removal in single-cell RNA-seq analysis""",
    paper_url="https://www.nature.com/articles/s41467-020-15851-3",
    paper_year=2020,
    code_url="https://github.com/eleozzr/desc",
    code_version=check_version("desc"),
    image="openproblems-python37-scgen"  # only if required
    # image="openproblems-python-batch-integration" # only if required
)
def desc_hvg_unscaled(adata):
    from ._utils import hvg_batch
    from scanpy.preprocessing import neighbors
    from scib.integration import runDESC

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata.obsm["X_emb"] = runDESC(adata, "batch").obsm["X_emb"]
    neighbors(adata, use_rep="X_emb")
    return adata


@method(
    method_name="DESC (hvg/scaled)",
    paper_name="""Deep learning enables accurate clustering with
     batch effect removal in single-cell RNA-seq analysis""",
    paper_url="https://www.nature.com/articles/s41467-020-15851-3",
    paper_year=2020,
    code_url="https://github.com/eleozzr/desc",
    code_version=check_version("desc"),
    image="openproblems-python37-scgen"  # only if required
    # image="openproblems-python-batch-integration" # only if required
)
def desc_hvg_scaled(adata):
    from ._utils import hvg_batch
    from scanpy.preprocessing import neighbors
    from scib.integration import runDESC
    from ._utils import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata.obsm["X_emb"] = runDESC(adata, "batch").obsm["X_emb"]
    neighbors(adata, use_rep="X_emb")
    return adata


@method(
    method_name="DESC (full/scaled)",
    paper_name="""Deep learning enables accurate clustering with
     batch effect removal in single-cell RNA-seq analysis""",
    paper_url="https://www.nature.com/articles/s41467-020-15851-3",
    paper_year=2020,
    code_url="https://github.com/eleozzr/desc",
    code_version=check_version("desc"),
    image="openproblems-python37-scgen"  # only if required
    # image="openproblems-python-batch-integration" # only if required
)
def desc_full_scaled(adata):
    from scanpy.preprocessing import neighbors
    from scib.integration import runDESC
    from ._utils import scale_batch

    adata = scale_batch(adata, "batch")
    adata.obsm["X_emb"] = runDESC(adata, "batch").obsm["X_emb"]
    neighbors(adata, use_rep="X_emb")
    return adata

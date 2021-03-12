# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version

import scprep

# from scIB.integration import _harmony


_harmony = scprep.run.RFunction(
    setup="""
            library(SingleCellExperiment)
            library(harmony)
            library(Seurat)
        """,
    args="sobj, batch",
    body="""
            seu <- as.Seurat(sobj)
            seu <- ScaleData(seu)
            seu <- RunPCA(seu, features=rownames(seu@assays$RNA))
            seu <- RunHarmony(seu, batch)
            seu <- as.SingleCellExperiment(seu)
            seu
        """,
)


@method(
    method_name="Harmony",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def harmony_full_unscaled(adata):
    from scIB.preprocessing import reduce_data

    adata = _harmony(adata, "batch")
    adata.obsm["X_emb"] = adata.obsm["HARMONY"]
    reduce_data(adata, use_rep="X_emb")
    # Complete the result in-place
    return adata


@method(
    method_name="Harmony (hvg/unscaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    image="openproblems-r-scib",  # only if required
)
def harmony_hvg_unscaled(adata):
    from _hvg import hvg_batch
    from scIB.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = _harmony(adata, "batch")
    adata.obsm["X_emb"] = adata.obsm["HARMONY"]
    reduce_data(adata, use_rep="X_emb")
    return adata


@method(
    method_name="Harmony (hvg/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    image="openproblems-r-scib",  # only if required
)
def harmony_hvg_scaled(adata):
    from _hvg import hvg_batch
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = _harmony(adata, "batch")
    adata.obsm["X_emb"] = adata.obsm["HARMONY"]
    reduce_data(adata, use_rep="X_emb")
    return adata


@method(
    method_name="Harmony(full/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    image="openproblems-r-scib",  # only if required
)
def harmony_full_scaled(adata):
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch

    adata = scale_batch(adata, "batch")
    adata = _harmony(adata, "batch")
    adata.obsm["X_emb"] = adata.obsm["HARMONY"]
    reduce_data(adata, use_rep="X_emb")
    return adata

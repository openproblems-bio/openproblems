# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version

import scprep

# from scIB.integration import _liger


_liger = scprep.run.RFunction(
    setup="""
            library(SingleCellExperiment)
            library(liger)
            library(Seurat)
        """,
    args="sobj_uni, batch",
    body="""
            # Only counts is converted to liger object. To pass our own normalized data,
            # store it in the "counts" slot
            sobj <- as.Seurat(sobj_uni)
            sobj@assays$RNA@counts = sobj@assays$RNA@data

            # Create Liger object
            lobj = seuratToLiger(sobj, combined.seurat=T, meta.var=batch, renormalize=F,
                         remove.missing=F)

            # We only pass nomarlized data, so store it as such
            lobj@norm.data <- lobj@raw.data

            # Assign hvgs
            lobj@var.genes <- rownames(sobj@assays$RNA)

            lobj <- scaleNotCenter(lobj, remove.missing=F)

            # Use tutorial coarse k suggests of 20.
            lobj <- optimizeALS(lobj, k=20, thresh=5e-5, nrep=3)

            lobj <- quantileAlignSNF(lobj, resolution=0.4,
                 small.clust.thresh=20)

            # Store embedding in initial Seurat object
            # Code taken from ligerToSeurat() function from LIGER
            inmf.obj <- new(
                Class = "DimReduc", feature.loadings = t(lobj@W),
                cell.embeddings = lobj@H.norm, key = "X_emb"
            )
            sobj@reductions["X_emb"] <- inmf.obj
            sobj <- as.SingleCellExperiment(sobj)
            return(sobj)
        """,
)


@method(
    method_name="Liger",
    paper_name="""Single-Cell Multi-omic Integration Compares and
                  Contrasts Features of Brain Cell Identity""",
    paper_url="https://www.cell.com/cell/fulltext/S0092-8674%2819%2930504-5",
    paper_year=2019,
    code_url="https://github.com/welch-lab/liger",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def liger_full_unscaled(adata):
    from scIB.preprocessing import reduce_data

    adata = _liger(adata, "batch")
    adata.obsm["X_emb"] = adata.obsm["X_EMB"]
    reduce_data(adata, umap=False, use_rep="X_emb")
    # Complete the result in-place
    return adata


@method(
    method_name="Liger (hvg/unscaled)",
    paper_name="""Single-Cell Multi-omic Integration Compares and
                  Contrasts Features of Brain Cell Identity""",
    paper_url="https://www.cell.com/cell/fulltext/S0092-8674%2819%2930504-5",
    paper_year=2019,
    code_url="https://github.com/welch-lab/liger",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def liger_hvg_unscaled(adata):
    from ._hvg import hvg_batch
    from scIB.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = _liger(adata, "batch")
    adata.obsm["X_emb"] = adata.obsm["X_EMB"]
    reduce_data(adata, umap=False, use_rep="X_emb")
    return adata


@method(
    method_name="Liger (hvg scaled)",
    paper_name="""Single-Cell Multi-omic Integration Compares and
                  Contrasts Features of Brain Cell Identity""",
    paper_url="https://www.cell.com/cell/fulltext/S0092-8674%2819%2930504-5",
    paper_year=2019,
    code_url="https://github.com/welch-lab/liger",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def liger_hvg_scaled(adata):
    from ._hvg import hvg_batch
    from ._hvg import scale_batch
    from scIB.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = _liger(adata, "batch")
    adata.obsm["X_emb"] = adata.obsm["X_EMB"]
    reduce_data(adata, umap=False, use_rep="X_emb")
    return adata


@method(
    method_name="Liger (full/scaled)",
    paper_name="""Single-Cell Multi-omic Integration Compares and
                  Contrasts Features of Brain Cell Identity""",
    paper_url="https://www.cell.com/cell/fulltext/S0092-8674%2819%2930504-5",
    paper_year=2019,
    code_url="https://github.com/welch-lab/liger",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def liger_full_scaled(adata):
    from ._hvg import scale_batch
    from scIB.preprocessing import reduce_data

    adata = scale_batch(adata, "batch")
    adata = _liger(adata, "batch")
    adata.obsm["X_emb"] = adata.obsm["X_EMB"]
    reduce_data(adata, umap=False, use_rep="X_emb")
    return adata

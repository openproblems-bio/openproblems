from .....tools.decorators import method
from .....tools.utils import check_version

import scprep

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
            sobj <- as.Seurat(sobj_uni, data=NULL)
            sobj@assays$RNA <- sobj@assays$originalexp
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

            # Return embedding
            return(lobj@H.norm)
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
def liger_full_unscaled(adata, test=False):
    from scanpy.pp import neighbors

    adata.obsm["X_emb"] = _liger(adata, "batch")
    neighbors(adata, use_rep="X_emb")

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
def liger_hvg_unscaled(adata, test=False):
    from ._utils import hvg_batch
    from scanpy.pp import neighbors

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata.obsm["X_emb"] = _liger(adata, "batch")
    neighbors(adata, use_rep="X_emb")
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
def liger_hvg_scaled(adata, test=False):
    from ._utils import hvg_batch
    from ._utils import scale_batch
    from scanpy.pp import neighbors

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata.obsm["X_emb"] = _liger(adata, "batch")
    neighbors(adata, use_rep="X_emb")
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
def liger_full_scaled(adata, test=False):
    from ._utils import scale_batch
    from scanpy.pp import neighbors

    adata = scale_batch(adata, "batch")
    adata.obsm["X_emb"] = _liger(adata, "batch")
    neighbors(adata, use_rep="X_emb")
    return adata

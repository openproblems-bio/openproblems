# from ....tools.normalize import log_cpm
from .....tools.decorators import method

from .....tools.utils import check_version
#from scIB.integration import _liger
from scIB.preprocessing import hvg_batch, scale_batch, reduce_data
import scprep

_liger = scprep.run.RFunction(
        setup="""
            library(SingleCellExperiment)
            library(liger)
            library(Seurat)
        """,
        args="sobj, batch",
        body="""
            # Only counts is converted to liger object. To pass our own normalized data,
            # store it in the "counts" slot
            sobj@assays$RNA@counts = sobj@assays$RNA@data

            # Create Liger object
            lobj = seuratToLiger(sobj, combined.seurat=T, meta.var=batch, renormalize=F,
                         remove.missing=F)

            # We only pass nomarlized data, so store it as such
            lobj@norm.data <- lobj@raw.data

            # Assign hvgs
            lobj@var.genes <- hvg

            lobj <- scaleNotCenter(lobj, remove.missing=F) # Can"t do our own scaling atm

            # Use tutorial coarse k suggests of 20.
            lobj <- optimizeALS(lobj, k=k, thresh=5e-5, nrep=3)

            lobj <- quantileAlignSNF(lobj, resolution=res, small.clust.thresh=small.clust.thresh)

            # Store embedding in initial Seurat object
            # Code taken from ligerToSeurat() function from LIGER
            inmf.obj <- new(
                Class = "DimReduc", feature.loadings = t(lobj@W),
                cell.embeddings = lobj@H.norm, key = "X_emb"
            )
            sobj@reductions["X_emb"] <- inmf.obj
            return(sobj)
        """
)


@method(
    method_name="Liger",
    paper_name="Sc"
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    # image="openproblems-template-image" # only if required
)
def liger_full_unscaled(adata):
    adata = _liger(adata, "batch")
    reduce_data(adata, use_emb="X_emb")
    # Complete the result in-place
    return adata

def liger_hvg_unscaled(adata):
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)    
    adata = _liger(adata, "batch")
    reduce_data(adata, use_emb="X_emb")
    return adata

def liger_hvg_scaled(adata):
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)    
    adata = scale_batch(adata, "batch")
    adata = _liger(adata, "batch")
    reduce_data(adata, use_emb="X_emb")
    return adata

def liger_full_scaled(adata):
    adata = scale_batch(adata, "batch")
    adata = _liger(adata, "batch")
    reduce_data(adata, use_emb="X_emb")
    return adata



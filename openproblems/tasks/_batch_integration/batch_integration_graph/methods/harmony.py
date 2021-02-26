# from ....tools.normalize import log_cpm
from .....tools.decorators import method

from .....tools.utils import check_version
#from scIB.integration import _harmony
from scIB.preprocessing import hvg_batch, scale_batch, reduce_data
import scprep

_harmony = scprep.run.RFunction(
        setup="""
            library(SingleCellExperiment)
            library(harmony)
            library(Seurat)
        """,
        args="sobj, batch",
        body="""
            sobj <- ScaleData(sobj)
	        sobj <- RunPCA(sobj, features=rownames(sobj@assays$RNA))
	        sobj <- RunHarmony(sobj, batch)
	        sobj[['X_emb']] <- sobj[['harmony']]

           return(sobj)
        """
)


@method(
    method_name="FastMNN feature",
    paper_name="Sc"
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    # image="openproblems-template-image" # only if required
)
def harmony_full_unscaled(adata):
    adata = _harmony(adata, "batch")
    reduce_data(adata, use_emb='X_emb')
    # Complete the result in-place
    return adata

def harmony_hvg_unscaled(adata):
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)    
    adata = _harmony(adata, "batch")
    reduce_data(adata, use_emb='X_emb')
    return adata

def harmony_hvg_scaled(adata):
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)    
    adata = scale_batch(adata, 'batch')
    adata = _harmony(adata, "batch")
    reduce_data(adata, use_emb='X_emb')
    return adata

def harmony_full_scaled(adata):
    adata = scale_batch(adata, 'batch')
    adata = _harmony(adata, "batch")
    reduce_data(adata, use_emb='X_emb')
    return adata



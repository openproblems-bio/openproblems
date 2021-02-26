# from ....tools.normalize import log_cpm
from .....tools.decorators import method

from .....tools.utils import check_version
#from scIB.integration import _seurat
from scIB.preprocessing import hvg_batch, scale_batch, reduce_data
import scprep

_seurat = scprep.run.RFunction(
        setup="""
            library(SingleCellExperiment)
            library(Seurat)
        """,
        args="sce, batch, hvg=2000"
        body="""
            batch_list = SplitObject(data, split.by = batch)

            batch_list <- lapply(X = batch_list, FUN = function(x) {
		        x  <- ScaleData(x, features = hvg)
		        x <- RunPCA(x, features = hvg)
		        return(x)
	         })

	        anchors = FindIntegrationAnchors(
	            object.list = batch_list,
	            anchor.features = hvg,
 		        scale = T,
		        l2.norm = T,
		        dims = 1:30,
        	    k.anchor = 5,
        	    k.filter = 200,
        	    k.score = 30,
		        reduction = "rpca",
        	    max.features = 200,
        	    eps = 0)
	        integrated = IntegrateData(
        	    anchorset = anchors,
		        new.assay.name = "integrated",
        	    features = NULL,
        	    features.to.integrate = NULL,
        	    dims = 1:30,
        	    k.weight = 100,
        	    weight.reduction = NULL,
        	    sd.weight = 1,
        	    sample.tree = NULL,
        	    preserve.order = F,
        	    do.cpp = T,
        	    eps = 0,
        	    verbose = T)
	        return(integrated)
      """
)


@method(
    method_name="Seurat RPCA",
    paper_name="Sc"
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("seurat"),
    # image="openproblems-template-image" # only if required
)
def seuratrpca_full_unscaled(adata):
    adata = _seurat(adata, "batch")
    reduce_data(adata)
    # Complete the result in-place
    return adata

def seuratrpca_hvg_unscaled(adata):
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)    
    adata = _seurat(adata, "batch")
    reduce_data(adata)
    return adata

def seuratrpca_hvg_scaled(adata):
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)    
    adata = scale_batch(adata, 'batch')
    adata = _seurat(adata, "batch")
    reduce_data(adata)
    return adata

def seuratrpca_full_scaled(adata):
    adata = scale_batch(adata, 'batch')
    adata = _seurat(adata, "batch")
    reduce_data(adata)
    return adata



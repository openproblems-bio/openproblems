# from ....tools.normalize import log_cpm
from .....tools.decorators import method

from .....tools.utils import check_version
#from scIB.integration import _fastmnn_feature
from scIB.preprocessing import hvg_batch, scale_batch, reduce_data
import scprep

_fastmnn_feature = scprep.run.RFunction(
        setup="""
            library(SingleCellExperiment)
            library(batchelor)
        """,
        args="sobj, batch, hvg=2000",
        body="""
            expr <- sobj@assays$RNA@data

	        sce <- fastMNN(expr, batch = sobj@meta.data[[batch]])

	        sobj@assays$RNA <- CreateAssayObject(assay(sce, "reconstructed"))
	        sobj[['X_emb']] <- CreateDimReducObject(reducedDim(sce, "corrected"), key='fastmnn_')

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
def fastmnn_feature_full_unscaled(adata):
    adata = _fastmnn_feature(adata, "batch")
    reduce_data(adata)
    # Complete the result in-place
    return adata

def fastmnn_feature_hvg_unscaled(adata):
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)    
    adata = _fastmnn_feature(adata, "batch")
    reduce_data(adata)
    return adata

def fastmnn_feature_hvg_scaled(adata):
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)    
    adata = scale_batch(adata, 'batch')
    adata = _fastmnn_feature(adata, "batch")
    reduce_data(adata)
    return adata

def fastmnn_feature_full_scaled(adata):
    adata = scale_batch(adata, 'batch')
    adata = _fastmnn_feature(adata, "batch")
    reduce_data(adata)
    return adata



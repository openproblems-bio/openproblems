# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version

#from scIB.integration import _fastmnn_embed

import scprep

_fastmnn_embed = scprep.run.RFunction(
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
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    # image="openproblems-template-image" # only if required
)
def fastmnn_embed_full_unscaled(adata):
    from scIB.preprocessing import hvg_batch
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch
    adata = _fastmnn_embed(adata, "batch")
    reduce_data(adata, use_emb='X_emb')
    # Complete the result in-place
    return adata

@method(
    method_name="FastMNN feature (hvg/unscaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    # image="openproblems-template-image" # only if required
)
def fastmnn_embed_hvg_unscaled(adata):
    from scIB.preprocessing import hvg_batch
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = _fastmnn_embed(adata, "batch")
    reduce_data(adata, use_emb='X_emb')
    return adata

@method(
    method_name="FastMNN feature (hvg/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    # image="openproblems-template-image" # only if required
)
def fastmnn_embed_hvg_scaled(adata):
    from scIB.preprocessing import hvg_batch
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, 'batch')
    adata = _fastmnn_embed(adata, "batch")
    reduce_data(adata, use_emb='X_emb')
    return adata

@method(
    method_name="FastMNN feature (full/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    # image="openproblems-template-image" # only if required
)
def fastmnn_embed_full_scaled(adata):
    from scIB.preprocessing import hvg_batch
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch
    adata = scale_batch(adata, 'batch')
    adata = _fastmnn_embed(adata, "batch")
    reduce_data(adata, use_emb='X_emb')
    return adata

# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version

import scprep

# from scIB.integration import _seurat


_seurat = scprep.run.RFunction(
    setup="""
            library(SingleCellExperiment)
            library(Seurat)
        """,
    args="sce, batch, hvg=2000",
    body="""
            data <- as.Seurat(sce)
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
               eps = 0,
               verbose = T)
            as.SingleCellExperiment(integrated)
      """,
)


@method(
    method_name="Seurat RPCA",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def seuratrpca_full_unscaled(adata):
    from scIB.preprocessing import reduce_data

    adata = _seurat(adata, "batch")
    reduce_data(adata)
    # Complete the result in-place
    return adata


@method(
    method_name="Seurat RPCA (hvg/unscaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def seuratrpca_hvg_unscaled(adata):
    from scIB.preprocessing import hvg_batch
    from scIB.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = _seurat(adata, "batch")
    reduce_data(adata)
    return adata


@method(
    method_name="Seurat RPCA (hvg/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def seuratrpca_hvg_scaled(adata):
    from scIB.preprocessing import hvg_batch
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = _seurat(adata, "batch")
    reduce_data(adata)
    return adata


@method(
    method_name="Seurat RPCA (full/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def seuratrpca_full_scaled(adata):
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch

    adata = scale_batch(adata, "batch")
    adata = _seurat(adata, "batch")
    reduce_data(adata)
    return adata

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

            anchors = FindIntegrationAnchors(
                object.list = batch_list,
                anchor.features = hvg,
                scale = T,
                l2.norm = T,
                dims = 1:30,
                k.anchor = 5,
                k.filter = 200,
                k.score = 30,
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
    method_name="Seurat CCA",
    paper_name="Comprehensive Integration of Single-Cell Data",
    paper_url="https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8",
    paper_year=2019,
    code_url="https://satijalab.org/seurat/",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def seurat_full_unscaled(adata):
    from scIB.preprocessing import reduce_data

    adata = _seurat(adata, "batch")
    reduce_data(adata, umap=False)
    # Complete the result in-place
    return adata


@method(
    method_name="Seurat CCA (hvg/unscaled)",
    paper_name="Comprehensive Integration of Single-Cell Data",
    paper_url="https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8",
    paper_year=2019,
    code_url="https://satijalab.org/seurat/",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def seurat_hvg_unscaled(adata):
    from ._hvg import hvg_batch
    from scIB.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = _seurat(adata, "batch")
    reduce_data(adata, umap=False)
    return adata


@method(
    method_name="Seurat CCA (hvg/scaled)",
    paper_name="Comprehensive Integration of Single-Cell Data",
    paper_url="https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8",
    paper_year=2019,
    code_url="https://satijalab.org/seurat/",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def seurat_hvg_scaled(adata):
    from ._hvg import hvg_batch
    from ._hvg import scale_batch
    from scIB.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = _seurat(adata, "batch")
    reduce_data(adata, umap=False)
    return adata


@method(
    method_name="Seurat CCA (full/scaled)",
    paper_name="Comprehensive Integration of Single-Cell Data",
    paper_url="https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8",
    paper_year=2019,
    code_url="https://satijalab.org/seurat/",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def seurat_full_scaled(adata):
    from ._hvg import scale_batch
    from scIB.preprocessing import reduce_data

    adata = scale_batch(adata, "batch")
    adata = _seurat(adata, "batch")
    reduce_data(adata, umap=False)
    return adata

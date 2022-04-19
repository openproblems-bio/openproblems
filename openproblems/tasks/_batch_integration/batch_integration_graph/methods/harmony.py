# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version

import scprep

# from scib.integration import _harmony


_harmony = scprep.run.RFunction(
    setup="""
            library(SingleCellExperiment)
            library(harmony)
            library(Seurat)
        """,
    args="sobj, batch",
    body="""
            seu <- as.Seurat(sobj, data=NULL)
            seu <- ScaleData(seu)
            seu <- RunPCA(seu, features=rownames(seu@assays$originalexp))
            seu <- RunHarmony(seu, batch)
            seu <- as.SingleCellExperiment(seu)
            seu
        """,
)


@method(
    method_name="Harmony",
    paper_name="""Fast, sensitive and accurate integration
                  of single-cell data with Harmony""",
    paper_url="https://www.nature.com/articles/s41592-019-0619-0",
    paper_year=2019,
    code_url="https://github.com/immunogenomics/harmony",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def harmony_full_unscaled(adata):
    from scanpy.pp import neighbors

    adata = _harmony(adata, "batch")
    adata.obsm["X_emb"] = adata.obsm["HARMONY"]
    neighbors(adata, use_rep="X_emb")
    adata.obs.nFeature_originalexp = adata.obs.nFeature_originalexp.astype("int")
    # Complete the result in-place
    return adata


@method(
    method_name="Harmony (hvg/unscaled)",
    paper_name="""Fast, sensitive and accurate integration
                  of single-cell data with Harmony""",
    paper_url="https://www.nature.com/articles/s41592-019-0619-0",
    paper_year=2019,
    code_url="https://github.com/immunogenomics/harmony",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def harmony_hvg_unscaled(adata):
    from ._utils import hvg_batch
    from scanpy.pp import neighbors

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = _harmony(adata, "batch")
    adata.obsm["X_emb"] = adata.obsm["HARMONY"]
    adata.obs.nFeature_originalexp = adata.obs.nFeature_originalexp.astype("int")
    neighbors(adata, use_rep="X_emb")
    return adata


@method(
    method_name="Harmony (hvg/scaled)",
    paper_name="""Fast, sensitive and accurate integration
                  of single-cell data with Harmony""",
    paper_url="https://www.nature.com/articles/s41592-019-0619-0",
    paper_year=2019,
    code_url="https://github.com/immunogenomics/harmony",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def harmony_hvg_scaled(adata):
    from ._utils import hvg_batch
    from ._utils import scale_batch
    from scanpy.pp import neighbors

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = _harmony(adata, "batch")
    adata.obsm["X_emb"] = adata.obsm["HARMONY"]
    adata.obs.nFeature_originalexp = adata.obs.nFeature_originalexp.astype("int")
    neighbors(adata, use_rep="X_emb")
    return adata


@method(
    method_name="Harmony(full/scaled)",
    paper_name="""Fast, sensitive and accurate integration
                  of single-cell data with Harmony""",
    paper_url="https://www.nature.com/articles/s41592-019-0619-0",
    paper_year=2019,
    code_url="https://github.com/immunogenomics/harmony",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",  # only if required
)
def harmony_full_scaled(adata):
    from ._utils import scale_batch
    from scanpy.pp import neighbors

    adata = scale_batch(adata, "batch")
    adata = _harmony(adata, "batch")
    adata.obsm["X_emb"] = adata.obsm["HARMONY"]
    adata.obs.nFeature_originalexp = adata.obs.nFeature_originalexp.astype("int")
    neighbors(adata, use_rep="X_emb")
    return adata

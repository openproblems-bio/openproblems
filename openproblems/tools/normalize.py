from . import decorators

import anndata as ad
import logging
import scanpy as sc
import scprep

log = logging.getLogger("openproblems")

_scran = scprep.run.RFunction(
    setup="""
        library('scran')
        library('BiocParallel')
        """,
    args="sce, min.mean=0.1",
    body="""
    if(class(sce@assays@data$X)!="dgCMatrix"){ # dgRMatrix conversion bug
    sce@assays@data$X <- as(as.matrix(sce@assays@data$X), Class="sparseMatrix")
    }
    sce <- computeSumFactors(
           sce, min.mean=min.mean,
           assay.type="X",
           BPPARAM=SerialParam()
           )
    sizeFactors(sce)
    """,
)


@decorators.normalizer
def log_scran_pooling(adata: ad.AnnData) -> ad.AnnData:
    """Normalize data with scran via rpy2."""
    scprep.run.install_bioconductor("scran")
    adata.obs["size_factors"] = _scran(adata)
    adata.X = scprep.utils.matrix_vector_elementwise_multiply(
        adata.X, adata.obs["size_factors"], axis=0
    )
    sc.pp.log1p(adata)
    return adata


def _cpm(adata: ad.AnnData):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e6, key_added="size_factors")


@decorators.normalizer
def cpm(adata: ad.AnnData) -> ad.AnnData:
    """Normalize data to counts per million."""
    _cpm(adata)
    return adata


@decorators.normalizer
def log_cpm(adata: ad.AnnData) -> ad.AnnData:
    """Normalize data to log counts per million."""
    _cpm(adata)
    sc.pp.log1p(adata)
    return adata


@decorators.normalizer
def sqrt_cpm(adata: ad.AnnData) -> ad.AnnData:
    """Normalize data to sqrt counts per million."""
    _cpm(adata)
    adata.X = scprep.transform.sqrt(adata.X)
    return adata


@decorators.normalizer
def log_cpm_hvg(adata: ad.AnnData, n_genes: int = 1000) -> ad.AnnData:
    """Normalize logCPM HVG

    Normalize data to log counts per million and select n_genes highly
    variable genes
    """

    adata = log_cpm(adata)

    if adata.n_vars < n_genes:
        log.warning(
            f"Less than {n_genes} genes, setting 'n_genes' to {int(adata.n_vars * 0.5)}"
        )
        n_genes = int(adata.n_vars * 0.5)

    sc.pp.highly_variable_genes(adata, n_top_genes=n_genes, flavor="cell_ranger")
    adata = adata[:, adata.var["highly_variable"]].copy()

    return adata

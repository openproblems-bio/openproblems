import anndata as ad
import scanpy as sc
import scprep


_scran = scprep.run.RFunction(
    setup="""
        library('scran')
        library('BiocParallel')
        """,
    args="sce, min.mean=0.1",
    body="""
    sce <- computeSumFactors(
           sce, min.mean=min.mean,
           assay.type="X",
           BPPARAM=BiocParallel::MulticoreParam()
           )
    sizeFactors(sce)
    """
)


def log_scran_pooling(adata: ad.AnnData) -> ad.AnnData:
    """Normalize data with scran via rpy2."""
    adata.obs["size_factors"] = _scran(adata)
    adata.X = scprep.utils.matrix_vector_elementwise_multiply(
        adata.X, adata.obs["size_factors"], axis=0
    )
    sc.pp.log1p(adata)
    return adata


def _cpm(adata: ad.AnnData):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e6, key_added="size_factors")


def log_cpm(adata: ad.AnnData) -> ad.AnnData:
    """Normalize data to log counts per million."""
    _cpm(adata)
    sc.pp.log1p(adata)
    return adata

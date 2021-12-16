import scanpy as sc
import scprep


def log_scran_pooling(adata):
    """Normalize data with scran via rpy2."""
    _scran = scprep.run.RFunction(
        setup="library('scran')",
        args="sce, min.mean=0.1",
        body="""
        sce <- computeSumFactors(
            sce, min.mean=min.mean,
            assay.type="X"
        )
        sizeFactors(sce)
        """,
    )
    adata.obs["size_factors"] = _scran(adata)
    adata.X = scprep.utils.matrix_vector_elementwise_multiply(
        adata.X, adata.obs["size_factors"], axis=0
    )
    sc.pp.log1p(adata)

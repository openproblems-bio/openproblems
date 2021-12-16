import scanpy as sc
import scprep
import scib


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


def scale_batch(adata, batch_key):
    """
    Scale count matrix by batch
    """
    return scib.pp.scale_batch(adata, batch=batch_key)


# TODO: HVG by batch
def hvg_batch(adata, batch_key, n_hvg):
    """
    Compute highly variable genes by batch
    """
    if n_hvg > adata.n_vars:
        return adata.var_names.tolist()
    return scib.pp.hvg_batch(
        adata,
        batch_key=batch_key,
        target_genes=n_hvg,
        adataOut=False
    )

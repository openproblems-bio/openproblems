from . import decorators
from anndata import AnnData

import scanpy as sc
import scprep

_scran = scprep.run.RFunction(
    setup="library('scran')",
    args="sce, min.mean=0.1",
    body="""
    sce <- computeSumFactors(sce, min.mean=min.mean, assay.type="X")
    sizeFactors(sce)
    """,
)


@decorators.normalizer
def log_scran_pooling(adata):
    """Normalize data with scran via rpy2."""
    scprep.run.install_bioconductor("scran")
    adata.obs["size_factors"] = _scran(adata)
    adata.X = scprep.utils.matrix_vector_elementwise_multiply(
        adata.X, adata.obs["size_factors"], axis=0
    )
    sc.pp.log1p(adata)


def _cpm(adata):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e6, key_added="size_factors")


@decorators.normalizer
def cpm(adata):
    """Normalize data to counts per million."""
    _cpm(adata)


@decorators.normalizer
def log_cpm(adata):
    """Normalize data to log counts per million."""
    _cpm(adata)
    sc.pp.log1p(adata)


@decorators.normalizer
def sqrt_cpm(adata):
    """Normalize data to sqrt counts per million."""
    _cpm(adata)
    adata.X = scprep.transform.sqrt(adata.X)


@decorators.normalizer
def preprocess_logCPM_1kHVG(adata: AnnData) -> AnnData:

    log_cpm(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=1000, flavor="cell_ranger")
    adata = adata[:, adata.var["highly_variable"]]
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    adata.obsm["X_input"] = adata.obsm["X_pca"]

    return adata

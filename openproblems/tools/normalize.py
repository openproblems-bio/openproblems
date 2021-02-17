from . import decorators
from . import utils

import scanpy as sc
import scprep


def _scran(adata, precluster=True):
    import anndata2ri
    import scIB.preprocessing

    # deactivate converter
    anndata2ri.deactivate()

    scIB.preprocessing.normalize(adata)

    try:
        utils.assert_finite(adata.X)
    except AssertionError:
        if precluster:
            adata.X = adata.raw.X
            _scran(adata, precluster=False)
        else:
            raise

    # deactivate converter
    anndata2ri.deactivate()

    # Make lightweight
    del adata.raw


@decorators.normalizer
def log_scran_pooling(adata):
    """Normalize data with scran via rpy2."""
    scprep.run.install_bioconductor("scran")
    _scran(adata)


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

from . import decorators
from . import utils
import numpy as np
import scanpy as sc
import scprep
import warnings


def _scran(adata, retries=2):
    import anndata2ri
    import scIB.preprocessing

    try:
        # Normalize via scran-pooling with own clustering at res=0.5
        if np.any(adata.X.sum(axis=1) < 1):
            warnings.warn("adata has empty cells")
        scIB.preprocessing.normalize(adata)
        utils.assert_finite(adata.X)
    except Exception as e:
        if np.all(adata.X.sum(axis=1) > 0):
            warnings.warn("adata has no empty cells")
        if retries > 0:
            warnings.warn(
                "scran pooling failed with {}({})".format(type(e).__name__, str(e))
            )
            _scran(adata, retries=retries - 1)
        else:
            raise e

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

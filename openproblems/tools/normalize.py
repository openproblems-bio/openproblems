import scanpy as sc
import scprep
import scIB.preprocessing
from .utils import normalizer


@normalizer
def log_scran_pooling(adata):
    """
    This function scran-normalizes the data
    """
    scprep.run.install_bioconductor("scran")
    # Normalize via scran-pooling with own clustering at res=0.5
    scIB.preprocessing.normalize(adata)

    # Make lightweight
    del adata.raw


def _cpm(adata):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e6, key_added="size_factors")


@normalizer
def cpm(adata):
    _cpm(adata)


@normalizer
def log_cpm(adata):
    _cpm(adata)
    sc.pp.log1p(adata)

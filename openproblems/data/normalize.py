import scanpy as sc
from scIB.preprocessing import normalize


def log_scran_pooling(adata):
    """
    This function scran-normalizes the data 
    """

    # Normalize via scran-pooling with own clustering at res=0.5
    normalize(adata)

    # Make lightweight
    del adata.raw


def cpm(adata):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e6, key_added="size_factors")


def log_cpm(adata):
    cpm(adata)
    sc.pp.log1p(adata)

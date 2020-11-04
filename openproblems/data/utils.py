import scanpy as sc
import os
import anndata
import hashlib

from decorator import decorator
from . import TEMPDIR


def _func_to_bytes(func):
    return bytes(".".join([func.__module__, func.__name__]), encoding="utf-8")


def _obj_to_bytes(obj):
    return bytes(str(obj), encoding="utf-8")


@decorator
def loader(func, *args, **kwargs):
    hash = hashlib.sha256()
    hash.update(_func_to_bytes(func))
    hash.update(_obj_to_bytes(args))
    hash.update(_obj_to_bytes(kwargs))
    filename = "openproblems_{}.h5ad".format(hash.hexdigest())
    filepath = os.path.join(TEMPDIR, filename)
    if os.path.isfile(filepath):
        return anndata.read_h5ad(filepath)
    else:
        adata = func(*args, **kwargs)
        adata.write(filepath)
        return adata


def filter_genes_cells(adata):
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_genes=1)


def subsample_even(adata, n_obs, even_obs):
    """Subsample a dataset evenly across an obs.

    Parameters
    ----------
    adata : AnnData
    n_obs : int
        Total number of cells to retain
    even_obs : str
        `adata.obs[even_obs]` to be subsampled evenly across partitions.

    Returns
    -------
    adata : AnnData
        Subsampled AnnData object
    """
    obs = adata.obs[even_obs].unique()
    adata_sub = []
    for o in obs:
        adata_sub.append(
            sc.pp.subsample(adata[adata.obs[even_obs] == o], n_obs=n_obs // len(obs))
        )
    return sc.concat(adata_sub)

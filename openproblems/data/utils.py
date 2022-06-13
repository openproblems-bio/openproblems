from . import TEMPDIR

import anndata
import functools
import hashlib
import logging
import os
import scanpy as sc

log = logging.getLogger("openproblems")


def _func_to_bytes(func):
    return bytes(".".join([func.__module__, func.__name__]), encoding="utf-8")


def _obj_to_bytes(obj):
    return bytes(str(obj), encoding="utf-8")


def _hash_function(func, *args, **kwargs):
    hash = hashlib.sha256()
    hash.update(_func_to_bytes(func))
    hash.update(_obj_to_bytes(args))
    hash.update(_obj_to_bytes(kwargs))
    return hash.hexdigest()


def _cache_path(func, *args, **kwargs):
    if hasattr(func, "__wrapped__"):
        func = func.__wrapped__
    filename = "openproblems_{}.h5ad".format(_hash_function(func, *args, **kwargs))
    return os.path.join(TEMPDIR, filename)


def loader(data_url):
    """Decorate a data loader function.

    Parameters
    ----------
    data_url : str
        Link to the original source of the dataset
    """

    def decorator(func):
        @functools.wraps(func)
        def apply_func(*args, **kwargs):
            filepath = _cache_path(func, *args, **kwargs)
            if os.path.isfile(filepath):
                log.debug(
                    "Loading cached {}({}, {}) dataset".format(
                        func.__name__, args, kwargs
                    )
                )
                adata = anndata.read_h5ad(filepath)
                adata.uns["_from_cache"] = True
                return adata
            else:
                log.debug(
                    "Downloading {}({}, {}) dataset".format(func.__name__, args, kwargs)
                )
                adata = func(*args, **kwargs)
                adata.uns["_from_cache"] = False
                if "counts" not in adata.layers:
                    adata.layers["counts"] = adata.X
                try:
                    os.mkdir(TEMPDIR)
                except OSError:
                    pass
                adata.write_h5ad(filepath)
                return adata

        apply_func.metadata = dict(data_url=data_url)
        return apply_func

    return decorator


def filter_genes_cells(adata):
    """Remove empty cells and genes."""
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_counts=2)


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
    values = adata.obs[even_obs].unique()
    adata_out = None
    n_obs_per_value = n_obs // len(values)
    for v in values:
        adata_subset = adata[adata.obs[even_obs] == v].copy()
        sc.pp.subsample(adata_subset, n_obs=min(n_obs_per_value, adata_subset.shape[0]))
        if adata_out is None:
            adata_out = adata_subset
        else:
            adata_out = adata_out.concatenate(adata_subset, batch_key="_obs_batch")

    adata_out.uns = adata.uns
    adata_out.varm = adata.varm
    adata_out.varp = adata.varp
    return adata_out

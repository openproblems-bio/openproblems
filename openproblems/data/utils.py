import scanpy as sc
import os
import anndata
import hashlib
import functools

from decorator import decorator
from . import TEMPDIR


def _func_to_bytes(func):
    return bytes(".".join([func.__module__, func.__name__]), encoding="utf-8")


def _obj_to_bytes(obj):
    return bytes(str(obj), encoding="utf-8")


def _bytes_to_str(series):
    if series.dtype.name == "object":
        is_bytes = series.apply(lambda x: isinstance(x, bytes))
        series[is_bytes] = series[is_bytes].astype(str)
    elif series.dtype.name == "category":
        series = series.astype(str).astype("category")
    return series


def _adata_bytes_to_str(adata):
    for meta in [adata.obs, adata.var]:
        meta.index = _bytes_to_str(meta.index.to_series())
        for key in meta.columns:
            meta[key] = _bytes_to_str(meta[key])
    for key, value in adata.uns.items():
        if isinstance(value, pd.Series):
            adata.uns[key] = _bytes_to_str(value)
    return adata


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
        adata = _adata_bytes_to_str(adata)
        try:
            os.mkdir(TEMPDIR)
        except OSError:
            pass
        adata.write(filepath)
        return adata


def filter_genes_cells(adata):
    """Remove empty cells and genes."""
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

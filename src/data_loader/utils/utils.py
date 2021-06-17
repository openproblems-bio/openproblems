import tempfile
import hashlib
import logging
import os
import scanpy as sc

log = logging.getLogger("openproblems")


def _make_tempdir():
    tempdir = os.path.join(tempfile.gettempdir(), "openproblems_cache")
    try:
        os.mkdir(tempdir)
        log.debug("Created data cache directory")
    except OSError:
        log.debug("Data cache directory exists")
    return tempdir


TEMPDIR = _make_tempdir()


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


def filter_genes_cells(adata):
    """Remove empty cells and genes."""
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_counts=2)

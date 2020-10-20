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

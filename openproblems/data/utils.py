import scanpy as sc
import os
import anndata

from decorator import decorator
from . import TEMPDIR


@decorator
def loader(func, *args, **kwargs):
    filename = "openproblems_{}-{}-{}.h5ad".format(
        hash(func), hash(str(args)), hash(str(kwargs))
    )
    filepath = os.path.join(TEMPDIR.name, filename)
    if os.path.isfile(filepath):
        return anndata.read_h5ad(filepath)
    else:
        adata = func(*args, **kwargs)
        adata.write(filepath)
        return adata


def filter_genes_cells(adata):
    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.filter_cells(adata, min_counts=1)

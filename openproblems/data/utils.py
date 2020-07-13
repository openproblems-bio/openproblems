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

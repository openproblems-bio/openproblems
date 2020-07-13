import os
import anndata
from decorator import decorator
from . import TEMPDIR


@decorator
def loader(func, *args, test=False, **kwargs):
    filename = "openproblems_{}_test-{}.h5ad".format(hash(func), test)
    filepath = os.path.join(TEMPDIR.name, filename)
    if os.path.isfile(filepath):
        return anndata.read_h5ad(filename)
    else:
        adata = func(test=test)
        adata.write(filename)
        return adata

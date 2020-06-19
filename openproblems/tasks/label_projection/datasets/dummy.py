import numpy as np
import anndata


def dummy():
    return anndata.AnnData(np.random.uniform((100, 10)))

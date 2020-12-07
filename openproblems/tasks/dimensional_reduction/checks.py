import numpy as np


def check_method(adata):

    assert np.sum(adata.obsm["X_emb"]) > 0

    return True

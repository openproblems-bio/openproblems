import numpy as np


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert np.sum(adata.obsm["X_emb"]) > 0
    return True

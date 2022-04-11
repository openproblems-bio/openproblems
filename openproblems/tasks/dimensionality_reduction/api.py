from ...data.sample import load_sample_data

import numpy as np
import scanpy as sc


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "X_emb" in adata.obsm
    assert adata.obsm["X_emb"].shape[1] == 2
    assert np.all(np.isfinite(adata.obsm["X_emb"]))
    return True


def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    return load_sample_data()


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    sc.tl.pca(adata)
    adata.obsm["X_emb"] = adata.obsm["X_pca"][:, :2]
    return adata

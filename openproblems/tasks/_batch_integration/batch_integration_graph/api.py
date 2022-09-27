from ....data.sample import load_sample_data
from ....tools.decorators import dataset

import numpy as np
import scanpy as sc


def check_dataset(adata):
    """Check that dataset output fits expected API."""

    assert "X_uni" in adata.obsm
    assert "batch" in adata.obs
    assert "labels" in adata.obs
    assert "uni_connectivities" in adata.obsp
    assert "log_normalized" in adata.layers

    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "connectivities" in adata.obsp
    assert "distances" in adata.obsp
    return True


@dataset()
def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    adata = load_sample_data()
    adata.obsm["X_uni"] = sc.pp.pca(adata.X)
    adata.obs["batch"] = np.random.choice(2, adata.shape[0], replace=True).astype(str)
    adata.obs["labels"] = np.random.choice(5, adata.shape[0], replace=True).astype(str)

    sc.pp.neighbors(adata, use_rep="X_uni", key_added="uni")
    return adata


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    import scanpy as sc

    sc.pp.neighbors(adata, use_rep="X_uni")
    return adata

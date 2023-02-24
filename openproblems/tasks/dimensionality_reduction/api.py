from ...data.sample import load_sample_data
from ...tools.decorators import dataset
from ...tools.normalize import log_cp10k

import numpy as np


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    assert "n_genes" in adata.uns
    assert adata.uns["n_genes"] == adata.shape[1]
    assert "X_ranking" in adata.obsm
    assert adata.obsm["X_ranking"].shape == (adata.shape[0], adata.shape[0])
    return True


def check_method(adata, is_baseline=False):
    """Check that method output fits expected API."""
    # check adata.X has not changed
    assert adata.uns["n_genes"] == adata.shape[1]
    assert adata.X is adata.layers["log_cp10k"]
    # check output
    assert "X_emb" in adata.obsm
    if not is_baseline:
        assert adata.obsm["X_emb"].shape[1] == 2
    assert np.all(np.isfinite(adata.obsm["X_emb"]))
    return True


@dataset()
def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    adata = load_sample_data()
    adata = log_cp10k(adata)
    adata.uns["n_genes"] = adata.shape[1]
    return adata


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    import scanpy as sc

    sc.tl.pca(adata)
    adata.obsm["X_emb"] = adata.obsm["X_pca"][:, :2]
    return adata

from ...data.sample import load_sample_data

import numpy as np


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    assert "train" in adata.obsm
    assert "test" in adata.obsm
    assert adata.obsm["train"].shape == adata.X.shape
    assert adata.obsm["test"].shape == adata.X.shape
    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "denoised" in adata.obsm
    assert adata.obsm["denoised"].shape == adata.X.shape
    return True


def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    adata = load_sample_data()
    adata.obsm["train"] = adata.X.toarray()
    adata.obsm["train"] = np.random.binomial(
        n=adata.obsm["train"].astype(int), p=0.8, size=adata.obsm["train"].shape
    ).astype(float)
    adata.obsm["test"] = adata.X.toarray()
    adata.obsm["test"].data = np.random.binomial(
        n=adata.obsm["test"].astype(int), p=0.2, size=adata.obsm["test"].shape
    ).astype(float)
    return adata


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    adata.obsm["denoised"] = adata.X.toarray() * 0.2
    return adata

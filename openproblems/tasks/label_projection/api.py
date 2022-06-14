from ...data.sample import load_sample_data

import numpy as np


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    assert "labels" in adata.obs
    assert "batch" in adata.obs
    assert "is_train" in adata.obs
    assert np.sum(adata.obs["is_train"]) > 0
    assert np.sum(~adata.obs["is_train"]) > 0
    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "labels_pred" in adata.obs
    return True


def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    adata = load_sample_data()
    adata.obs["batch"] = np.random.choice(2, adata.shape[0], replace=True).astype(str)
    adata.obs["labels"] = np.random.choice(5, adata.shape[0], replace=True).astype(str)
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True, p=[0.8, 0.2]
    )
    return adata


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    adata.obs["labels_pred"] = adata.obs["labels"]
    adata.obs.loc[adata.obs.index[::5], "labels_pred"] = np.random.choice(
        5, len(adata.obs["labels_pred"][::5]), replace=True
    ).astype(str)
    return adata

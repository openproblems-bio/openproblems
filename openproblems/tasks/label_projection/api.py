from ...data.sample import load_sample_data
from ...tools.decorators import dataset

import numpy as np
import pandas as pd


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    assert "labels" in adata.obs
    assert "batch" in adata.obs
    assert "is_train" in adata.obs
    assert np.issubdtype(adata.obs["is_train"].dtype, bool)
    assert pd.api.types.is_categorical(adata.obs["batch"])
    assert pd.api.types.is_categorical(adata.obs["labels"])
    assert np.sum(adata.obs["is_train"]) > 0
    assert np.sum(~adata.obs["is_train"]) > 0
    return True


def check_method(adata, is_baseline=False):
    """Check that method output fits expected API."""
    assert "labels_pred" in adata.obs
    return True


@dataset()
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
    adata.obs["labels_pred"][::5] = np.random.choice(
        adata.obs["labels"].cat.categories,
        len(adata.obs["labels_pred"][::5]),
        replace=True,
    )
    return adata

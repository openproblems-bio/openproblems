from ....data.sample import load_sample_data
from ....tools.decorators import dataset

import numpy as np
import scanpy as sc


def check_dataset(adata):
    """Check that dataset output fits expected API."""

    assert "X_uni_pca" in adata.obsm
    assert "batch" in adata.obs
    assert "labels" in adata.obs
    assert "log_normalized" in adata.layers

    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "X_emb" in adata.obsm
    return True


@dataset()
def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    adata = load_sample_data()

    adata.var.index = adata.var.gene_short_name.astype(str)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata["log_normalized"] = adata.X
    adata.obsm["X_uni_pca"] = sc.pp.pca(adata.X)
    adata.obs["batch"] = np.random.choice(2, adata.shape[0], replace=True).astype(str)
    adata.obs["labels"] = np.random.choice(5, adata.shape[0], replace=True).astype(str)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""

    adata.obsm["X_emb"] = adata.obsm["X_uni_pca"]
    return adata

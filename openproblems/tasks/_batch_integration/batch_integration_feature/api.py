from ....data.sample import load_sample_data
from ....tools.decorators import dataset
import numpy as np
import scanpy as sc


def check_dataset(adata):
    """Check that dataset output fits expected API."""

    assert "batch" in adata.obs
    assert "labels" in adata.obs
    assert "log_normalized" in adata.layers
    assert "counts" in adata.layers
    assert adata.var_names.is_unique
    assert adata.obs_names.is_unique

    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "log_normalized" in adata.layers
    assert adata.layers["log_normalized"] is not adata.X
    return True


@dataset()
def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    adata = load_sample_data()

    adata.var.index = adata.var.gene_short_name.astype(str)
    sc.pp.normalize_total(adata)

    adata.obsm["X_uni"] = sc.pp.pca(adata.X)
    adata.obs["batch"] = np.random.choice(2, adata.shape[0], replace=True).astype(str)
    adata.obs["labels"] = np.random.choice(5, adata.shape[0], replace=True).astype(str)
    adata.layers["counts"] = adata.X
    adata.layers["log_normalized"] = adata.X.multiply(
        10000 / adata.X.sum(axis=1)
    ).tocsr()
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    adata.X = adata.X.multiply(2)
    return adata

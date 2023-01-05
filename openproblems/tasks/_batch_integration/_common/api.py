from ....data.sample import load_sample_data
from ....tools.decorators import dataset
from .utils import filter_celltypes

import numpy as np

MIN_CELLS_PER_CELLTYPE = 50


def check_dataset(adata):
    """Check that dataset output fits expected API."""

    assert "batch" in adata.obs
    assert "labels" in adata.obs
    assert (adata.obs["labels"].value_counts() >= MIN_CELLS_PER_CELLTYPE).all()

    assert "log_normalized" in adata.layers
    return True


@dataset()
def sample_dataset(run_pca: bool = False, run_neighbors: bool = False):
    """Create a simple dataset to use for testing methods in this task."""
    import scanpy as sc

    adata = load_sample_data()
    adata.uns["organism"] = "human"

    adata.var.index = adata.var.gene_short_name.astype(str)
    sc.pp.normalize_total(adata)
    adata.layers["log_normalized"] = adata.X

    adata.obs["batch"] = np.random.choice(2, adata.shape[0], replace=True).astype(str)
    adata.obs["labels"] = np.random.choice(5, adata.shape[0], replace=True).astype(str)
    adata = filter_celltypes(adata)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    if run_pca:
        adata.obsm["X_uni_pca"] = sc.pp.pca(adata.X)
    if run_neighbors:
        sc.pp.neighbors(adata, use_rep="X_uni_pca", key_added="uni")
    return adata

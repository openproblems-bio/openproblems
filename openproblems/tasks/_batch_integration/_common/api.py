from ....data.sample import load_sample_data
from ....tools.decorators import dataset
from .utils import filter_celltypes
from .utils import precompute_hvg

import numbers
import numpy as np

MIN_CELLS_PER_CELLTYPE = 50
N_HVG_UNINT = 2000


def check_neighbors(adata, neighbors_key, connectivities_key, distances_key):
    assert neighbors_key in adata.uns
    assert adata.uns[neighbors_key]["connectivities_key"] == connectivities_key
    assert adata.uns[neighbors_key]["distances_key"] == distances_key
    assert connectivities_key in adata.obsp
    assert distances_key in adata.obsp


def check_dataset(
    adata,
    do_check_hvg=False,
):
    """Check that dataset output fits expected API."""

    assert "batch" in adata.obs
    assert "labels" in adata.obs
    assert (adata.obs["labels"].value_counts() >= MIN_CELLS_PER_CELLTYPE).all()

    assert "log_normalized" in adata.layers
    assert "counts" in adata.layers

    assert adata.var_names.is_unique
    assert adata.obs_names.is_unique

    assert "n_genes_pre" in adata.uns
    assert isinstance(adata.uns["n_genes_pre"], numbers.Integral)
    assert adata.uns["n_genes_pre"] == adata.n_vars

    assert "organism" in adata.uns
    assert adata.uns["organism"] in ["mouse", "human"]

    assert "X_uni_pca" in adata.obsm

    if do_check_hvg:
        assert "hvg_unint" in adata.uns
        assert len(adata.uns["hvg_unint"]) == min(N_HVG_UNINT, adata.n_vars)
        assert np.all(np.isin(adata.uns["hvg_unint"], adata.var.index))

    check_neighbors(adata, "uni", "uni_connectivities", "uni_distances")

    return True


@dataset()
def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    import scanpy as sc

    adata = load_sample_data()
    adata.uns["organism"] = "human"

    adata.var.index = adata.var.gene_short_name.astype(str)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    sc.pp.normalize_total(adata)
    adata.layers["log_normalized"] = adata.X

    adata.obs["batch"] = np.random.choice(2, adata.shape[0], replace=True).astype(str)
    adata.obs["labels"] = np.random.choice(3, adata.shape[0], replace=True).astype(str)
    adata = filter_celltypes(adata)

    adata.uns["hvg_unint"] = precompute_hvg(adata)
    adata.uns["n_genes_pre"] = adata.n_vars

    adata.obsm["X_uni_pca"] = sc.pp.pca(adata.X)
    sc.pp.neighbors(adata, use_rep="X_uni_pca", key_added="uni")
    return adata

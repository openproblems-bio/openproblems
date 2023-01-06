from ....tools.decorators import dataset
from .._common import api

MIN_CELLS_PER_CELLTYPE = 50


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    api.check_dataset(adata)

    assert "X_uni_pca" in adata.obsm

    assert "uni" in adata.uns
    assert adata.uns["uni"]["connectivities_key"] == "uni_connectivities"
    assert adata.uns["uni"]["distances_key"] == "uni_distances"
    assert "uni_connectivities" in adata.obsp
    assert "uni_distances" in adata.obsp

    return True


def check_method(adata, is_baseline=False):
    """Check that method output fits expected API."""
    assert "neighbors" in adata.uns
    assert adata.uns["neighbors"]["connectivities_key"] == "connectivities"
    assert adata.uns["neighbors"]["distances_key"] == "distances"
    assert "connectivities" in adata.obsp
    assert "distances" in adata.obsp
    return True


@dataset()
def sample_dataset():
    return api.sample_dataset(run_pca=True, run_neighbors=True)


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    import scanpy as sc

    sc.pp.neighbors(adata, use_rep="X_uni_pca")
    return adata

from .datasets.immune import immune_batch


def check_dataset(adata):
    """Check that dataset output fits expected API."""

    assert "X_uni" in adata.obsm
    assert "batch" in adata.obs
    assert "labels" in adata.obs
    assert "uni_connectivities" in adata.obsp

    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "connectivities" in adata.obsp
    assert "distances" in adata.obsp
    return True


def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    adata = immune_batch(True)
    # print(adata.obs.columns)

    return adata


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    import scanpy as sc

    sc.pp.neighbors(adata, use_rep="X_uni")
    return adata

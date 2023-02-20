from .._common import api

MIN_CELLS_PER_CELLTYPE = 50

check_dataset = api.check_dataset
sample_dataset = api.sample_dataset


def check_method(adata, is_baseline=False):
    """Check that method output fits expected API."""
    api.check_neighbors(adata, "neighbors", "connectivities", "distances")
    return True


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    import scanpy as sc

    sc.pp.neighbors(adata, use_rep="X_uni_pca")
    return adata

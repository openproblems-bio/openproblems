from ....tools.decorators import dataset
from .._common import api

import functools

MIN_CELLS_PER_CELLTYPE = 50


check_dataset = functools.partial(
    api.check_dataset, do_check_pca=True, do_check_neighbors=True
)


def check_method(adata, is_baseline=False):
    """Check that method output fits expected API."""
    api.check_neighbors(adata, "neighbors", "connectivities", "distances")
    return True


@dataset()
def sample_dataset():
    return api.sample_dataset(run_pca=True, run_neighbors=True)


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    import scanpy as sc

    sc.pp.neighbors(adata, use_rep="X_uni_pca")
    return adata

from ....tools.decorators import dataset
from .._common import api

import functools

check_dataset = functools.partial(api.check_dataset, do_check_pca=True)


def check_method(adata, is_baseline=False):
    """Check that method output fits expected API."""
    assert "X_emb" in adata.obsm
    # check organism was not removed
    assert "organism" in adata.uns
    return True


@dataset()
def sample_dataset():
    return api.sample_dataset(run_pca=True)


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""

    adata.obsm["X_emb"] = adata.obsm["X_uni_pca"]
    return adata

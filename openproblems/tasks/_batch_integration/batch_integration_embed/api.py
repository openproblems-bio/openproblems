from .._common import api

check_dataset = api.check_dataset
sample_dataset = api.sample_dataset


def check_method(adata, is_baseline=False):
    """Check that method output fits expected API."""
    assert "X_emb" in adata.obsm
    # check organism was not removed
    assert "organism" in adata.uns
    return True


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""

    adata.obsm["X_emb"] = adata.obsm["X_uni_pca"]
    adata.uns["is_baseline"] = False
    return adata

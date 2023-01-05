from .._common import api


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    api.check_dataset(adata)

    assert "counts" in adata.layers
    assert adata.var_names.is_unique
    assert adata.obs_names.is_unique

    return True


def check_method(adata, is_baseline=False):
    """Check that method output fits expected API."""
    assert "log_normalized" in adata.layers
    if not is_baseline:
        assert adata.layers["log_normalized"] is not adata.X
    return True


sample_dataset = api.sample_dataset


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    adata.X = adata.X.multiply(2)
    return adata

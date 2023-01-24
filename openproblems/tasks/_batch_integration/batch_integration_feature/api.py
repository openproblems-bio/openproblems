from .._common import api

check_dataset = api.check_dataset


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

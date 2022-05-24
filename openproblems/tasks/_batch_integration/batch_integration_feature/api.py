from ..batch_integration_graph.datasets.pancreas import pancreas_batch


def check_dataset(adata):
    """Check that dataset output fits expected API."""

    assert "X_uni" in adata.obsm
    assert "batch" in adata.obs
    assert "labels" in adata.obs

    return True


def check_method(adata):
    """Check that method output fits expected API."""
    return True


def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    adata = pancreas_batch(True)

    return adata


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""

    return adata

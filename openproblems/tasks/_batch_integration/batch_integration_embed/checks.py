def check_dataset(adata):
    """Check that dataset output fits expected API."""
    assert "batch" in adata.obs
    assert "labels" in adata.obs
    return True


def check_method(adata):
    """Check that method output fits expected API."""
    # assert "template_output" in adata.obs
    assert "distances" in adata.obsp
    return True

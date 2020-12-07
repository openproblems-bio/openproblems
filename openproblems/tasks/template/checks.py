def check_dataset(adata):
    """Check that dataset output fits expected API."""
    # TODO: update
    assert "template_variable" in adata.obs
    return True


def check_method(adata):
    """Check that method output fits expected API."""
    # TODO: update
    assert "template_output" in adata.obs
    return True

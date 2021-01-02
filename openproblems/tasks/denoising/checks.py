def check_dataset(adata):
    """Check that dataset output fits expected API."""
    assert "train" in adata.obsm
    assert "test" in adata.obsm
    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "denoised" in adata.obsm
    return True

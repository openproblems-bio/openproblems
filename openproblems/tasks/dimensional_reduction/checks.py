def check_dataset(adata):
    """Check that dataset output fits expected API."""
    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "X_emb" in adata.obsm
    assert adata.obsm["X_emb"].shape[1] == 2
    return True

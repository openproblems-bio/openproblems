import numpy as np


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    for key in ["species", "release"]:
        assert key in adata.uns, key
        assert isinstance(adata.uns[key], str), key

    assert "mode2" in adata.obsm
    assert "mode2_obs" in adata.uns
    assert np.all(adata.obs.index == adata.uns["mode2_obs"])

    for key in [
        "mode2_var",
        "mode2_var_chr",
        "mode2_var_start",
        "mode2_var_end",
    ]:
        assert key in adata.uns, key
        assert len(adata.uns[key]) == adata.obsm["mode2"].shape[1], key

    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "gene_score" in adata.obsm
    return True

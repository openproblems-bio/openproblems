from ...data.multimodal.sample import load_sample_data
from ...tools.decorators import dataset

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
    assert adata.obsm["gene_score"].shape == adata.X.shape
    return True


@dataset()
def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    adata = load_sample_data()

    adata.uns["species"] = "mus_musculus"
    adata.uns["version"] = "GRCm38"
    adata.uns["release"] = "100"

    return adata


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    adata.obsm["gene_score"] = adata.X.toarray() / adata.X.max() + np.random.normal(
        0, 0.1, adata.X.shape
    )
    return adata

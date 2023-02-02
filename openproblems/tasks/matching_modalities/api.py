from ...data.multimodal.sample import load_sample_data
from ...tools.decorators import dataset

import numpy as np


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    assert "mode2" in adata.obsm
    assert "mode2_obs" in adata.uns
    assert "mode2_var" in adata.uns
    assert np.all(adata.obs.index == adata.uns["mode2_obs"])
    assert len(adata.uns["mode2_var"]) == adata.obsm["mode2"].shape[1]
    return True


def check_method(adata, is_baseline=False):
    """Check that method output fits expected API."""
    assert "aligned" in adata.obsm
    assert "mode2_aligned" in adata.obsm
    assert adata.obsm["aligned"].shape[0] == adata.shape[0]
    assert adata.obsm["mode2_aligned"].shape[0] == adata.obsm["mode2"].shape[0]
    assert adata.obsm["aligned"].shape[1] == adata.obsm["mode2_aligned"].shape[1]
    assert np.all(np.isfinite(adata.obsm["aligned"]))
    assert np.all(np.isfinite(adata.obsm["mode2_aligned"]))
    return True


@dataset()
def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    return load_sample_data()


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    adata.obsm["aligned"] = adata.X / adata.X.max() + np.random.normal(
        0, 0.1, adata.X.shape
    )
    adata.obsm["mode2_aligned"] = adata.X / adata.X.max() + np.random.normal(
        0, 0.1, adata.X.shape
    )
    return adata

import numpy as np


def check_dataset(adata):
    assert "mode2" in adata.obsm
    assert "mode2_raw" in adata.obsm
    assert "mode2_obs" in adata.uns
    assert "mode2_var" in adata.uns
    assert np.all(adata.obs.index == adata.uns["mode2_obs"])
    assert len(adata.uns["mode2_var"]) == adata.obsm["mode2"].shape[1]
    return True


def check_method(adata):
    assert "aligned" in adata.obsm
    assert "mode2_aligned" in adata.obsm
    assert adata.obsm["aligned"].shape[0] == adata.shape[0]
    assert adata.obsm["mode2_aligned"].shape[0] == adata.obsm["mode2"].shape[0]
    assert adata.obsm["aligned"].shape[1] == adata.obsm["mode2_aligned"].shape[1]
    return True

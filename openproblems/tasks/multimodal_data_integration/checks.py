import numpy as np


def check_dataset(adata):
    assert "mode2" in adata.uns
    assert adata.shape[0] == adata.uns["mode2"].shape[0]
    assert np.all(adata.obs.index == adata.uns["mode2"].obs.index)
    assert adata.raw is not None
    assert adata.uns["mode2"].raw is not None
    return True


def check_method(adata):
    assert "aligned" in adata.obsm
    assert "aligned" in adata.uns["mode2"].obsm
    assert adata.obsm["aligned"].shape[0] == adata.shape[0]
    assert adata.uns["mode2"].obsm["aligned"].shape[0] == adata.uns["mode2"].shape[0]
    assert adata.obsm["aligned"].shape[1] == adata.uns["mode2"].obsm["aligned"].shape[1]
    return True

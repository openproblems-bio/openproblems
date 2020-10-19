import numpy as np


def check_dataset(adata):
    assert "mode2" in adata.obsm
    assert "mode2_obs" in adata.uns
    assert "mode2_var" in adata.uns
    assert np.all(adata.obs.index == adata.uns["mode2_obs"])
    assert len(adata.uns["mode2_var"]) == adata.obsm["mode2"].shape[1]
    return True


def check_method(adata):
    assert "atac_rna_cor" in adata.obs.columns
    return True

import numpy as np


def check_dataset(adata):
    assert "labels" in adata.obs
    assert "is_train" in adata.obs
    assert np.sum(adata.obs["is_train"]) > 0
    assert np.sum(~adata.obs["is_train"]) > 0
    return True


def check_method(adata):
    assert "labels_pred" in adata.obs
    return True

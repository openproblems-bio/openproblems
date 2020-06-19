def check_dataset(adata):
    assert "labels" in adata.obs
    assert "is_train" in adata.obs
    return True


def check_method(adata):
    assert "labels_pred" in adata.obs
    return True

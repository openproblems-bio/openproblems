import numpy as np


def check_dataset(adata):
    assert "train" in adata.obsm
    assert "test" in adata.obsm
    return True


def check_method(adata):
    assert "denoised" in adata.obsm
    return True

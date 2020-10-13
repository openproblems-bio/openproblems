import numpy as np


def check_dataset(adata):
    assert "template_variable" in adata.obs
    return True


def check_method(adata):
    assert "template_output" in adata.obs
    return True

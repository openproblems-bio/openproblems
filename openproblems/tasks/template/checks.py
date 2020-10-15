import numpy as np


def check_dataset(adata):
    # TODO: update
    assert "template_variable" in adata.obs
    return True


def check_method(adata):
    # TODO: update
    assert "template_output" in adata.obs
    return True

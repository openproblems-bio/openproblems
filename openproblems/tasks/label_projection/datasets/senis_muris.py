import numpy as np
from ....data.senis_muris import load_senis_muris


def senis_muris_timepoint(test=False):
    adata = load_senis_muris(test=test)
    adata.obs["is_train"] = adata.obs["timepoint"]
    return adata


def senis_muris_batch(test=False):
    adata = load_senis_muris(test=test)
    adata.obs["is_train"] = adata.obs["batch"]
    return adata


def senis_muris_random(test=False):
    adata = load_senis_muris(test=test)
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True
    )
    return adata

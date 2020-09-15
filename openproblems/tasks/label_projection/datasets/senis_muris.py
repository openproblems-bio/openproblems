import numpy as np
from ....data.senis_muris import load_senis_muris
from ....tools.decorators import dataset


@dataset("Tabula Senis (by timepoint)")
def senis_muris_timepoint(test=False):
    adata = load_senis_muris(test=test)
    adata.obs["is_train"] = adata.obs["age"] != "3m"
    return adata


@dataset("Tabula Senis (random split)")
def senis_muris_random(test=False):
    adata = load_senis_muris(test=test)
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True
    )
    return adata

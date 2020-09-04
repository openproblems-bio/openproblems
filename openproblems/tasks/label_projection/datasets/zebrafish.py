import numpy as np
from ....data.zebrafish import load_zebrafish


def zebrafish_labels(test=False):
    adata = load_zebrafish(test=test)
    adata.obs["labels"] = adata.obs["cell_type"]
    adata.obs["is_train"] = adata.obs["lab"] == "Schier"
    return adata


def zebrafish_random(test=False):
    adata = load_zebrafish(test=test)
    adata.obs["labels"] = adata.obs["cell_type"]
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True
    )
    return adata

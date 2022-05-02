from ....data.zebrafish import load_zebrafish
from ....tools.decorators import dataset

import numpy as np


@dataset("Zebrafish (by labels)")
def zebrafish_labels(test=False):
    adata = load_zebrafish(test=test)
    adata.obs["labels"] = adata.obs["cell_type"]
    adata.obs["batch"] = adata.obs["lab"]
    adata.obs["is_train"] = adata.obs["lab"] == adata.obs["lab"][0]
    return adata


@dataset("Zebrafish (random split)")
def zebrafish_random(test=False):
    adata = load_zebrafish(test=test)
    adata.obs["labels"] = adata.obs["cell_type"]
    adata.obs["batch"] = adata.obs["lab"]
    adata.obs["is_train"] = np.random.choice(
        [True, False], adata.shape[0], replace=True, p=[0.8, 0.2]
    )
    return adata

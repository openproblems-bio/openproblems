import numpy as np


def dummy(adata):
    adata.obs["labels_pred"] = np.random.choice(2, adata.shape[0], replace=True)

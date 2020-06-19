import numpy as np


def accuracy(adata):
    return np.mean(adata.obs["labels"] == adata.obs["labels_pred"])

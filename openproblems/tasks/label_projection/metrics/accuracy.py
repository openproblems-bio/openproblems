import numpy as np


def accuracy(adata):
    test_data = adata[~adata.obs["is_train"]]
    return np.mean(test_data.obs["labels"] == test_data.obs["labels_pred"])

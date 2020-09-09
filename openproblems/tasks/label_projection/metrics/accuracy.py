import numpy as np

from ....tools.decorators import metric


@metric(metric_name="Accuracy", maximize=True)
def accuracy(adata):
    test_data = adata[~adata.obs["is_train"]]
    return np.mean(test_data.obs["labels"] == test_data.obs["labels_pred"])

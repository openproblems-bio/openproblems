import numpy as np
from sklearn.preprocessing import LabelEncoder

from ....tools.decorators import metric


@metric(metric_name="Accuracy", maximize=True)
def accuracy(adata):
    encoder = LabelEncoder().fit(adata.obs["labels"])
    test_data = adata[~adata.obs["is_train"]]

    test_data.obs["labels"] = encoder.transform(test_data.obs["labels"])
    test_data.obs["labels_pred"] = encoder.transform(test_data.obs["labels_pred"])

    return np.mean(test_data.obs["labels"] == test_data.obs["labels_pred"])

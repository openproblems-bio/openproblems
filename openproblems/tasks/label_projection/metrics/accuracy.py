from ....tools.decorators import metric

import numpy as np


@metric(metric_name="Accuracy", paper_reference="grandini2020metrics", maximize=True)
def accuracy(adata):
    import sklearn.preprocessing

    encoder = sklearn.preprocessing.LabelEncoder().fit(adata.obs["labels"])
    test_data = adata[~adata.obs["is_train"]]

    labels = encoder.transform(test_data.obs["labels"])
    labels_pred = encoder.transform(test_data.obs["labels_pred"])

    return np.mean(labels == labels_pred)

from sklearn.metrics import f1_score
from sklearn.preprocessing import LabelEncoder

from ....tools.decorators import metric


def _f1(adata, average="weighted"):
    encoder = LabelEncoder().fit(adata.obs["labels"])
    test_data = adata[~adata.obs["is_train"]]

    test_data.obs["labels"] = encoder.transform(test_data.obs["labels"])
    test_data.obs["labels_pred"] = encoder.transform(test_data.obs["labels_pred"])

    return f1_score(
        test_data.obs["labels"], test_data.obs["labels_pred"], average=average
    )


@metric(metric_name="F1 score", maximize=True)
def f1(adata):
    return _f1(adata, average="weighted")


@metric(metric_name="Micro F1 score", maximize=True)
def f1_micro(adata):
    return _f1(adata, average="micro")

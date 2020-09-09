from sklearn.metrics import f1_score

from ....tools.decorators import metric


@metric(metric_name="F1 score", maximize=True)
def f1(adata):
    test_data = adata[~adata.obs["is_train"]]
    return f1_score(
        test_data.obs["labels"], test_data.obs["labels_pred"], average="weighted"
    )

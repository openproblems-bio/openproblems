from ....tools.decorators import metric


def _f1(adata, average="weighted"):
    import sklearn.metrics
    import sklearn.preprocessing

    encoder = sklearn.preprocessing.LabelEncoder().fit(adata.obs["labels"])
    test_data = adata[~adata.obs["is_train"]]

    test_data.obs["labels"] = encoder.transform(test_data.obs["labels"])
    test_data.obs["labels_pred"] = encoder.transform(test_data.obs["labels_pred"])

    return sklearn.metrics.f1_score(
        test_data.obs["labels"], test_data.obs["labels_pred"], average=average
    )


@metric(
    metric_name="F1 score",
    metric_summary=(
        "The [F1"
        " score](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html)"  # noqa: E501
        " is a weighted average of the precision and recall over all class labels,"
        " where an F1 score reaches its best value at 1 and worst score at 0, where"
        " each class contributes to the score relative to its frequency in the dataset."
    ),
    paper_reference="grandini2020metrics",
    maximize=True,
)
def f1(adata):
    return _f1(adata, average="weighted")


@metric(
    metric_name="Macro F1 score",
    metric_summary=(
        "The macro F1 score is an unweighted F1 score, where each class "
        "contributes equally, regardless of its frequency."
    ),
    paper_reference="grandini2020metrics",
    maximize=True,
)
def f1_macro(adata):
    return _f1(adata, average="macro")

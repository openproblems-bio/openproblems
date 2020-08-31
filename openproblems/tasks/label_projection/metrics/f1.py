from sklearn.metrics import f1_score


def f1(adata):
    test_data = adata[~adata.obs["is_train"]]
    return f1_score(test_data.obs["labels"], test_data.obs["labels_pred"])

from sklearn.linear_model import LogisticRegression
import numpy as np
from ....tools.normalize import log_cpm, log_scran_pooling


def _logistic_regression(adata):
    classifier = LogisticRegression()

    adata_train = adata[adata.obs["is_train"]]
    adata_test = adata[~adata.obs["is_train"]].copy()

    classifier.fit(adata_train.X, adata_train.obs["labels"])
    adata_test.obs["labels_pred"] = classifier.predict(adata_test.X)

    adata.obs["labels_pred"] = [
        adata_test.obs["labels_pred"][idx] if idx in adata_test.obs_names else np.nan
        for idx in adata.obs_names
    ]


def logistic_regression_log_cpm(adata):
    log_cpm(adata)
    _logistic_regression(adata)


def logistic_regression_scran(adata):
    log_scran_pooling(adata)
    _logistic_regression(adata)

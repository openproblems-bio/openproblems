from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_scran_pooling
from ....tools.utils import check_version
from .utils import pca_op

import numpy as np
import scipy.sparse
import sklearn.linear_model
import sklearn.pipeline
import sklearn.preprocessing


def _logistic_regression(adata, max_iter=1000, n_pca=100):
    adata_train = adata[adata.obs["is_train"]]
    adata_test = adata[~adata.obs["is_train"]].copy()
    is_sparse = scipy.sparse.issparse(adata.X)

    classifier = sklearn.pipeline.Pipeline(
        [
            ("pca", pca_op(adata_train, adata_test, n_components=n_pca)),
            ("scaler", sklearn.preprocessing.StandardScaler(with_mean=not is_sparse)),
            ("regression", sklearn.linear_model.LogisticRegression(max_iter=max_iter)),
        ]
    )

    # Fit to train data
    classifier.fit(adata_train.X, adata_train.obs["labels"].astype(str))

    # Predict on test data
    adata_test.obs["labels_pred"] = classifier.predict(adata_test.X)

    adata.obs["labels_pred"] = [
        adata_test.obs["labels_pred"][idx] if idx in adata_test.obs_names else np.nan
        for idx in adata.obs_names
    ]
    return adata


@method(
    method_name="Logistic regression (log CPM)",
    paper_name="Applied Logistic Regression",
    paper_url="https://books.google.com/books?id=64JYAwAAQBAJ",
    paper_year=2013,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.linear_model.LogisticRegression.html",
    code_version=check_version("scikit-learn"),
)
def logistic_regression_log_cpm(adata):
    log_cpm(adata)
    return _logistic_regression(adata)


@method(
    method_name="Logistic regression (log scran)",
    paper_name="Applied Logistic Regression",
    paper_url="https://books.google.com/books?id=64JYAwAAQBAJ",
    paper_year=2013,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.linear_model.LogisticRegression.html",
    code_version=check_version("scikit-learn"),
    image="openproblems-r-base",
)
def logistic_regression_scran(adata):
    log_scran_pooling(adata)
    return _logistic_regression(adata)

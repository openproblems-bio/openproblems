from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

import numpy as np
from scipy import sparse

from ....tools.normalize import log_cpm, log_scran_pooling
from ....tools.decorators import method
from ....tools.utils import check_version


def _logistic_regression(adata, max_iter=1000, n_pca=100):

    adata_train = adata[adata.obs["is_train"]]
    adata_test = adata[~adata.obs["is_train"]].copy()
    is_sparse = sparse.issparse(adata.X)

    min_pca = min([adata_train.shape[0], adata_test.shape[0], adata.shape[1]])
    if is_sparse:
        min_pca -= 1
    n_pca = min([n_pca, min_pca])
    pca_op = TruncatedSVD if is_sparse else PCA

    classifier = Pipeline(
        [
            ("pca", pca_op(n_components=n_pca)),
            ("scaler", StandardScaler(with_mean=not is_sparse)),
            ("regression", LogisticRegression(max_iter=max_iter)),
        ]
    )

    # Fit to train data
    classifier.fit(adata_train.X, adata_train.obs["labels"])

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
    code_url="https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html",
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
    code_url="https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html",
    code_version=check_version("scikit-learn"),
    image="openproblems-r-base",
)
def logistic_regression_scran(adata):
    log_scran_pooling(adata)
    return _logistic_regression(adata)

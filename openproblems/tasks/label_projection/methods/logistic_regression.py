from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_scran_pooling
from ....tools.utils import check_version
from .sklearn import classifier

import sklearn.linear_model


def _logistic_regression(adata, test=False, max_iter=None):
    if test:
        max_iter = max_iter or 100
    else:  # pragma: no cover
        max_iter = max_iter or 1000
    return classifier(
        adata, estimator=sklearn.linear_model.LogisticRegression, max_iter=max_iter
    )


@method(
    method_name="Logistic regression (log CPM)",
    paper_name="Applied Logistic Regression",
    paper_url="https://books.google.com/books?id=64JYAwAAQBAJ",
    paper_year=2013,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.linear_model.LogisticRegression.html",
    code_version=check_version("scikit-learn"),
)
def logistic_regression_log_cpm(adata, test=False, max_iter=None):
    adata = log_cpm(adata)
    return _logistic_regression(adata, test=test, max_iter=max_iter)


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
def logistic_regression_scran(adata, test=False, max_iter=None):
    adata = log_scran_pooling(adata)
    return _logistic_regression(adata, test=test, max_iter=max_iter)

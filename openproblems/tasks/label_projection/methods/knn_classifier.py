from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_scran_pooling
from ....tools.utils import check_version

import sklearn.neighbors
from .sklearn import classifier


@method(
    method_name="K-neighbors classifier (log CPM)",
    paper_name="Nearest neighbor pattern classification",
    paper_url="https://doi.org/10.1109/TIT.1967.1053964",
    paper_year=1967,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.neighbors.KNeighborsClassifier.html",
    code_version=check_version("scikit-learn"),
)
def logistic_regression_log_cpm(adata):
    log_cpm(adata)
    return classifier(adata, estimator=sklearn.neighbors.KNeighborsClassifier)


@method(
    method_name="K-neighbors classifier (log scran)",
    paper_name="Nearest neighbor pattern classification",
    paper_url="https://doi.org/10.1109/TIT.1967.1053964",
    paper_year=1967,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.neighbors.KNeighborsClassifier.html",
    code_version=check_version("scikit-learn"),
    image="openproblems-r-base",
)
def logistic_regression_scran(adata):
    log_scran_pooling(adata)
    return classifier(adata, estimator=sklearn.neighbors.KNeighborsClassifier)

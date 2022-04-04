from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_scran_pooling
from ....tools.utils import check_version
from .sklearn import classifier

import sklearn.neighbors


@method(
    method_name="K-neighbors classifier (log CPM)",
    paper_name="Nearest neighbor pattern classification",
    paper_url="https://doi.org/10.1109/TIT.1967.1053964",
    paper_year=1967,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.neighbors.KNeighborsClassifier.html",
    code_version=check_version("scikit-learn"),
)
def knn_classifier_log_cpm(adata, test=False):
    adata = log_cpm(adata)
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
def knn_classifier_scran(adata, test=False):
    adata = log_scran_pooling(adata)
    return classifier(adata, estimator=sklearn.neighbors.KNeighborsClassifier)

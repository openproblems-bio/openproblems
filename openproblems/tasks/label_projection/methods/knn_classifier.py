from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_scran_pooling
from .sklearn import classifier

import functools

_knn_classifier_method = functools.partial(
    method,
    paper_name="Nearest neighbor pattern classification",
    paper_reference="cover1967nearest",
    paper_year=1967,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.neighbors.KNeighborsClassifier.html",
)


@_knn_classifier_method(
    method_name="K-neighbors classifier (log CPM)",
)
def knn_classifier_log_cpm(adata, test=False):
    import sklearn.neighbors

    adata = log_cpm(adata)
    return classifier(adata, estimator=sklearn.neighbors.KNeighborsClassifier)


@_knn_classifier_method(
    method_name="K-neighbors classifier (log scran)",
    image="openproblems-r-base",
)
def knn_classifier_scran(adata, test=False):
    import sklearn.neighbors

    adata = log_scran_pooling(adata)
    return classifier(adata, estimator=sklearn.neighbors.KNeighborsClassifier)

from ....tools.decorators import method
from ....tools.normalize import log_cp10k
from ....tools.normalize import log_scran_pooling
from .sklearn import classifier

import functools

_knn_classifier_method = functools.partial(
    method,
    method_summary=(
        'K-neighbors classifier uses the "k-nearest neighbours" approach, which is a'
        " popular machine learning algorithm for classification and regression tasks."
        " The assumption underlying KNN in this context is that cells with similar gene"
        " expression profiles tend to belong to the same cell type. For each unlabelled"
        " cell, this method computes the $k$ labelled cells (in this case, 5) with the"
        " smallest distance in PCA space, and assigns that cell the most common cell"
        " type among its $k$ nearest neighbors."
    ),
    paper_name="Nearest neighbor pattern classification",
    paper_reference="cover1967nearest",
    paper_year=1967,
    code_url=(
        "https://scikit-learn.org/stable/modules/generated/"
        "sklearn.neighbors.KNeighborsClassifier.html"
    ),
)


@_knn_classifier_method(
    method_name="K-neighbors classifier (log CP10k)",
)
def knn_classifier_log_cp10k(adata, test=False):
    import sklearn.neighbors

    adata = log_cp10k(adata)
    return classifier(adata, estimator=sklearn.neighbors.KNeighborsClassifier)


@_knn_classifier_method(
    method_name="K-neighbors classifier (log scran)",
    image="openproblems-r-base",
)
def knn_classifier_scran(adata, test=False):
    import sklearn.neighbors

    adata = log_scran_pooling(adata)
    return classifier(adata, estimator=sklearn.neighbors.KNeighborsClassifier)

from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_scran_pooling
from ....tools.utils import check_version

from .utils import pca_op

import numpy as np
import scipy.sparse

import sklearn.neural_network
import sklearn.pipeline
import sklearn.preprocessing


def _mlp(adata, n_pca=100):

    adata_train = adata[adata.obs["is_train"]]
    adata_test = adata[~adata.obs["is_train"]].copy()
    is_sparse = scipy.sparse.issparse(adata.X)

    classifier = sklearn.pipeline.Pipeline(
        [
            ("pca", pca_op(adata_train, adata_test, n_components=n_pca)),
            ("scaler", sklearn.preprocessing.StandardScaler(with_mean=not is_sparse)),
            ("regression", sklearn.neural_network.MLPClassifier()),
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
    method_name="Multilayer perceptron (log CPM)",
    paper_name="Connectionist learning procedures",
    paper_url="https://www.sciencedirect.com/science/article/pii/B9780080510552500298",
    paper_year=1990,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.neural_network.MLPClassifier.html",
    code_version=check_version("scikit-learn"),
)
def mlp_log_cpm(adata):
    log_cpm(adata)
    return _mlp(adata)


@method(
    method_name="Multilayer perceptron (log scran)",
    paper_name="Connectionist learning procedures",
    paper_url="https://www.sciencedirect.com/science/article/pii/B9780080510552500298",
    paper_year=1990,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.neural_network.MLPClassifier.html",
    code_version=check_version("scikit-learn"),
    image="openproblems-r-base",
)
def mlp_scran(adata):
    log_scran_pooling(adata)
    return _mlp(adata)

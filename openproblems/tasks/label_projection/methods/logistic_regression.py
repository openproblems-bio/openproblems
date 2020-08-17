from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
import sklearn
import numpy as np
from ....tools.normalize import log_cpm, log_scran_pooling


def _logistic_regression(adata):
    classifier = LogisticRegression(max_iter=1000)

    adata_train = adata[adata.obs["is_train"]]
    adata_test = adata[~adata.obs["is_train"]].copy()

    # Dimensionality reduction & scaling
    pc_op = PCA(n_components=100)
    adata_train.obsm['pca'] = pc_op.fit_transform(adata_train.X.toarray())
    adata_train.obsm['pca_scale'] = sklearn.preprocessing.scale(adata_train.obsm['pca'])

    # Fit to train data
    classifier.fit(adata_train.obsm['pca_scale'], adata_train.obs["labels"])

    # Reduce and scale test data
    adata_test.obsm['pca'] = adata_test.X @ pc_op.components_.T
    adata_test.obsm['pca_scale'] = sklearn.preprocessing.scale(adata_test.obsm['pca'])

    # Fit test data
    adata_test.obs["labels_pred"] = classifier.predict(adata_test.obsm['pca_scale'])

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

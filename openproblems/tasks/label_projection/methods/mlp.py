from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_scran_pooling
from ....tools.utils import check_version
from .sklearn import classifier

import sklearn.neural_network


def _mlp(adata, test=False):
    if test:
        hidden_layer_sizes = (20,)
        max_iter = 100
    else:
        hidden_layer_sizes = (100, 100)
        max_iter = 1000
    return classifier(
        adata,
        estimator=sklearn.neural_network.MLPClassifier,
        hidden_layer_sizes=hidden_layer_sizes,
        max_iter=max_iter
    )


@method(
    method_name="Multilayer perceptron (log CPM)",
    paper_name="Connectionist learning procedures",
    paper_url="https://www.sciencedirect.com/science/article/pii/B9780080510552500298",
    paper_year=1990,
    code_url="https://scikit-learn.org/stable/modules/generated/"
    "sklearn.neural_network.MLPClassifier.html",
    code_version=check_version("scikit-learn"),
)
def mlp_log_cpm(adata, test=False):
    log_cpm(adata)
    return _mlp(adata, test=test)


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
def mlp_scran(adata, test=False):
    log_scran_pooling(adata)
    return _mlp(adata, test=test)

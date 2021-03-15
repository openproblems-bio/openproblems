from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_scran_pooling
from ....tools.utils import check_version
from .sklearn import classifier

import sklearn.neural_network


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
    return classifier(adata, estimator=sklearn.neural_network.MLPClassifier)


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
    return classifier(adata, estimator=sklearn.neural_network.MLPClassifier)

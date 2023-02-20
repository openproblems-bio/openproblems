from ....tools.decorators import method
from ....tools.normalize import log_cp10k
from ....tools.normalize import log_scran_pooling
from .sklearn import classifier

import functools

_mlp_method = functools.partial(
    method,
    method_summary="TODO",
    paper_name="Connectionist learning procedures",
    paper_reference="hinton1989connectionist",
    paper_year=1990,
    code_url=(
        "https://scikit-learn.org/stable/modules/generated/"
        "sklearn.neural_network.MLPClassifier.html"
    ),
)


def _mlp(adata, test=False, max_iter=None, hidden_layer_sizes=None):
    import sklearn.neural_network

    if test:
        hidden_layer_sizes = hidden_layer_sizes or (20,)
        max_iter = max_iter or 100
    else:  # pragma: no cover
        hidden_layer_sizes = hidden_layer_sizes or (100, 100)
        max_iter = max_iter or 1000
    return classifier(
        adata,
        estimator=sklearn.neural_network.MLPClassifier,
        hidden_layer_sizes=hidden_layer_sizes,
        max_iter=max_iter,
    )


@_mlp_method(
    method_name="Multilayer perceptron (log CP10k)",
)
def mlp_log_cp10k(adata, test=False, max_iter=None, hidden_layer_sizes=None):
    adata = log_cp10k(adata)
    return _mlp(
        adata, test=test, max_iter=max_iter, hidden_layer_sizes=hidden_layer_sizes
    )


@_mlp_method(
    method_name="Multilayer perceptron (log scran)",
    image="openproblems-r-base",
)
def mlp_scran(adata, test=False, max_iter=None, hidden_layer_sizes=None):
    adata = log_scran_pooling(adata)
    return _mlp(
        adata, test=test, max_iter=max_iter, hidden_layer_sizes=hidden_layer_sizes
    )

from ....tools.decorators import method
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version
from anndata import AnnData
from scipy.sparse import issparse
from typing import Any
from typing import Dict
from typing import Tuple

import numpy as np

_CFG = {
    "hyperparameter": {
        "optimization": {"method": "Adam", "learning_rate": 0.01},
        "batch_size": 512,
        "max_epoch": 20,  # original: 100
        "regularizer_l2": 0.001,
        "perplexity": 10,
        "seed": 1,
        "normalize": None,
    },
    "architecture": {
        "latent_dimension": 2,
        "inference": {"layer_size": [128, 64, 32]},
        "model": {"layer_size": [32, 32, 32, 64, 128]},
        "activation": "ELU",
    },
}


def _init_model(
    x: np.ndarray, config: Dict[str, Any], is_train: bool = True
) -> Tuple[np.ndarray, None, Dict[str, Any], "DataSet", "SCVIS", float]:  # noqa: F821
    from scvis import data
    from scvis.model import SCVIS

    architecture = config["architecture"]
    architecture.update({"input_dimension": x.shape[1]})

    hyperparameter = config["hyperparameter"]
    if hyperparameter["batch_size"] > x.shape[0]:
        hyperparameter.update({"batch_size": x.shape[0]})

    model = SCVIS(architecture, hyperparameter)
    normalizer = 1.0

    if is_train:
        if hyperparameter["normalize"] is not None:
            normalizer = float(hyperparameter["normalize"])
        else:
            normalizer = np.max(np.abs(x))
    else:
        if hyperparameter["normalize"] is not None:
            normalizer = float(hyperparameter["normalize"])

    x /= normalizer

    y = None  # no labels
    train_data = data.DataSet(x, y)

    return x, y, hyperparameter, train_data, model, normalizer


def _fit(data: np.ndarray):
    x, y, hyperparameter, train_data, model, normalizer = _init_model(
        data, config=_CFG, is_train=True
    )

    iter_per_epoch = round(x.shape[0] / hyperparameter["batch_size"])
    max_iter = int(iter_per_epoch * hyperparameter["max_epoch"])
    # limit the max_iter because of CI
    max_iter = min(max_iter, 100)

    _ = model.train(
        data=train_data,
        batch_size=hyperparameter["batch_size"],
        verbose=False,
        verbose_interval=max_iter,
        show_plot=False,
        plot_dir=None,
        max_iter=max_iter,
        pretrained_model=None,
    )
    model.set_normalizer(normalizer)

    return model, x


@method(
    method_name="scvis (CPU) (logCPM, 1kHVG)",
    paper_name="Interpretable dimensionality reduction "
    "of single celltranscriptome data with deep generative models",
    paper_url="https://www.nature.com/articles/s41467-018-04368-5",
    paper_year=2018,
    code_url="https://bitbucket.org/jerry00/scvis-dev/",
    code_version=check_version("scvis"),
    image="openproblems-python36",
)
def scvis_logCM_1kHVG(adata: AnnData, test: bool = False) -> AnnData:
    adata = log_cpm_hvg(adata)

    model, x = _fit(adata.X.A if issparse(adata.X) else adata.X)
    emb, _ = model.encode(x)

    adata.obsm["X_emb"] = np.asarray(emb[:, :2])

    return adata
